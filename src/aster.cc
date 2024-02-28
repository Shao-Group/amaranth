/*
Part of Aster Transcript Assembler
(c) 2024 by Xiaofei Carl Zang, Mingfu Shao, and The Pennsylvania State University.
See LICENSE for licensing.
*/

#include <algorithm>
#include <cmath>
#include "util.h"
#include "aster.h"
#include "basic_algo.h"

aster::aster(const splice_graph &g, const hyper_set &h)
	: origr(g), gr(g), hs(h)
{
	paths.clear();					// predicted paths, original v index, inclusive
    trsts.clear();			// predicted transcripts, original v index
	non_full_trsts.clear();		// predicted non full length transcripts
	mode = asterMode;

    topological_sort_vertices();
	make_stats();
	init_edgeres();

	if(mode == aster_mode::STAT_ONLY) print_stats();

	if(verbose >= 1) cout << "aster assembling " << gr.gid << endl;
	assemble();
	get_transcripts();	
	assert(paths.size() == trsts.size() + non_full_trsts.size());
	if(verbose >= 1 && successStatus) cout << "aster assembled " << gr.gid << ", #path = " << paths.size() << endl;
}

int aster::assemble()
{	
	if(mode == aster_mode::STAT_ONLY) return 0;
	if(gr.num_edges() == 0) return 0;
	if(gr.num_vertices() == 2) return 0;
	assert(gr.num_vertices() > 2);

	if(gr.num_vertices() > 1000) //FIXME:
	{
		cerr << gr.gid <<" too big to process with D&C " << gr.num_vertices() << endl; 
		successStatus = false;
		return 0;
	}	
	if (true)	//CLEAN: balance?
	{
		for(int i = 1; i < gr.num_vertices() - 1; i++) balance_vertex(i);
		for(int i = 1; i < gr.num_vertices() - 1; i++) balance_vertex(i);
		for(int i = 1; i < gr.num_vertices() - 1; i++) balance_vertex(i);
		for(int i = 1; i < gr.num_vertices() - 1; i++) balance_vertex(i);
		remove_small_junctions();
		gr.refine_splice_graph();
	}

	try
	{
		divide_conquer();
		successStatus = true;
	}
	catch(const aster_error& e)
	{
		cerr << e.what() << '\n';
		print_stats();
		string gene_start_end = gr.chrm + ":"
							  + to_string(gr.get_vertex_info(0).lpos) + "-"
							  + to_string(gr.get_vertex_info(gr.num_vertices() - 1).rpos)
							  + "\n";
		gr.graphviz("asterviz_err." + gr.gid + ".dot", gene_start_end + tp2v_to_string());
		paths.clear();
		successStatus = false;
	}	
	
	if (verbose >= 2)  for(int i = 0; i < paths.size(); i++) paths[i].print(i);
	return 0;
}

int aster::divide_conquer()
{
	if(output_graphviz_files) 
	{
		string gene_start_end = gr.chrm + ":"
							  + to_string(gr.get_vertex_info(0).lpos) + "-"
							  + to_string(gr.get_vertex_info(gr.num_vertices() - 1).rpos)
							  + "\n";
		gr.graphviz("asterviz." + gr.gid + ".dot", gene_start_end + tp2v_to_string());
	}
	assert(gr.num_vertices() > 2);
	assert(tp2v.size() == gr.num_vertices());
	divide_conquer({tp2v});
	paths.clear();
	paths = res.subpaths;
	if(paths.size() >= 1) for(const path & p : paths)	assert(origr.valid_path(p.v));
	return 0;
}

// divide_conquer(i ,j) solves a subproblem between 
// i, j are tp2v indices
int aster::divide_conquer(aster_index ai)
{
	int s = ai.s();
	int t = ai.t();

	if(verbose >= 2)
	{
		string msg = "aster processing subgraph, vertex [" + to_string(s) + ", " + to_string(t) + "]"; 
		// msg += " (topoIndex [" + to_string(source) + "," + to_string(target) + "])";
		cout << msg << endl;
	}

	assert(s <= t);
	assert(s < gr.num_vertices() && t < gr.num_vertices() && s >= 0 && t >= 0);

	if (resolve_trivial_node(ai))		
	{
		dnc_counter_resolve_trivial_node ++;
		return 0;
	}
	if (divide_conquer_single_vertex(ai))			
	{
		dnc_counter_single ++;
		return 0;
	} 
	if (divide_conquer_unitig(ai))					
	{
		dnc_counter_unitig ++;
		return 0;
	}
	if (divide_conquer_abutting(ai))				
	{
		dnc_counter_abutting ++;
		return 0;
	}
	if (divide_conquer_cut_termini(ai)) 	
	{
		dnc_counter_cut_vertex ++;
		return 0;
	}
	if (divide_conquer_articulation_point(ai))		
	{
		dnc_counter_articulation_point_disjoint ++;
		return 0;
	}
	// if (resolve_trivial_intersection(ai))		
	// {
	// 	dnc_counter_resolve_trivial_intersection ++;
	// 	return 0;
	// }
	if (resolve_intersection_edge(ai))
	{
		dnc_counter_resolve_intersection_edge ++;
		return 0;
	}
	if(mode != aster_mode::MINI && greedy(ai))
	{
		return 0;
	}

	num_intersecting_graph ++;

	string msg = "aster-mini failed on graph " + gr.gid;
	msg += " [" + gr.chrm + ":" + to_string(gr.get_vertex_info(0).lpos) + "-" + to_string(gr.get_vertex_info(0).rpos) + "]";
	msg += " in subgraph vertexIndex [" + to_string(s) + ", " + to_string(t) + "]";
	throw aster_error(msg.c_str());
	return -1;
}

bool aster::greedy(aster_index ai)
{
	//TODO:
	return false;	
	throw runtime_error("not implemented yet"); 
}

bool aster::resolve_intersection_edge(aster_index ai)
{
	return false;
	divide_conquer(ai);
	return true;
}

bool aster::resolve_trivial_node(aster_index ai)
{
	int s = ai.s();
	int t = ai.t();
	assert(s < gr.num_vertices() && t < gr.num_vertices() && s >= 0 && t >= 0);
	assert(s < t);
	assert(gr.out_degree(s) >= 1 || gr.in_degree(t) >= 1);

	for(int i = 1; i < ai.size() - 1; i ++)
	{
		int v = ai.at(i);
		assert(v > s && v < t);
		if(gr.in_degree(v) > 1 && gr.out_degree(v) > 1) continue;

		PEEI inEdges = gr.in_edges(v);
		PEEI outEdges = gr.out_edges(v);
		edges_combine_consecutive_and_replace(inEdges, outEdges);
		assert(gr.degree(v) == 0);	

		if(verbose >= 2) 
		{
			string msg = "aster resolved a trivial node " + to_string(v); 
			msg += "in subgraph, vertex [" + to_string(s) + ", " + to_string(t) + "]"; 
			cout << msg << endl;
		}
		return true;
	}
	return false;
}

/* resolve trivila intersections only 
*  all four vertices should be directly connected but not recursively solvable
*  1->2->3->4 and 1->3 and 2->4 forming a diomond shape
*  s-> k1 and k2, k1 -> k2 and t, k2 -> t
*/
/* 
bool aster::resolve_trivial_intersection(int source, int target, aster_result& res) 
{
	assert(source < tp2v.size() && target < tp2v.size());
	int s = tp2v[source];
	int t = tp2v[target];
	assert(s < gr.num_vertices() && t < gr.num_vertices() && s >= 0 && t >= 0);
	assert(s < t - 2);
	assert(source < target - 2);
	assert(gr.out_degree(s) >= 1 || gr.in_degree(t) >= 1);
	
	if(gr.edge_exists(s, t)) return false;
	if(source + 3 != target) return false;

	int k1 = tp2v[source + 1];
	int k2 = tp2v[source + 2];
	assert(s < k1);
	assert(k1 < k2);
	assert(k2 < t);
	if (gr.out_degree(s) != 2) return false;
	if (gr.out_degree(k1)!= 2 || gr.in_degree(k1) != 1) return false;
	if (gr.out_degree(k2)!= 1 || gr.in_degree(k2) != 2) return false;
	if (gr.in_degree(t)  != 2) return false;
	if (!gr.edge_exists(s, k1)) return false;
	if (!gr.edge_exists(s, k2)) return false;
	if (!gr.edge_exists(k1, k2)) return false;
	if (!gr.edge_exists(k1, t)) return false;
	if (!gr.edge_exists(k2, t)) return false;

	double wsk1  = gr.get_edge_weight(gr.edge(s, k1).first);
	double wsk2  = gr.get_edge_weight(gr.edge(s, k2).first);
	double wk1k2 = gr.get_edge_weight(gr.edge(k1, k2).first);
	double wk1t  = gr.get_edge_weight(gr.edge(k1, t).first);
	double wk2t  = gr.get_edge_weight(gr.edge(k2, t).first);
	// average if true, geom mean if false
	bool   _avg_ = false;	
	double abd1 = _avg_? (wsk1 + wk1k2 + wk2t / 3.0) : pow(wsk1 * wk1k2 * wk2t, 1.0/3.0);
	double abd2 = _avg_? (wsk1 + wk1t / 2.0) : pow(wsk1 * wk1t, 1.0/2.0);
	double abd3 = _avg_? (wsk2 + wk2t / 2.0) : pow(wsk2 * wk2t, 1.0/2.0);
	res.subpaths.push_back({{s, k1, k2, t}, abd1});
	res.subpaths.push_back({{s, k1, t}, abd2});
	res.subpaths.push_back({{s, k2, t}, abd3});

	replace_aster_index_to_one_edge(source, target, abd1 + abd2 + abd3);

	if(verbose >= 2) 
	{
		string msg = "aster processed subgraph, vertex [" + to_string(s) + ", " + to_string(t) + "]"; 
		msg += " (topoIndex [" + to_string(source) + "," + to_string(target) + "]), ";
		msg += "a trivial intersecting graph"; 
		cout << msg << endl;
	}

	if(true) // CLEAN: debug trivial
	{
		string gene_start_end = gr.chrm + ":"
								+ to_string(gr.get_vertex_info(0).lpos) + "-"
								+ to_string(gr.get_vertex_info(gr.num_vertices() - 1).rpos)
								+ "\n";
		gr.graphviz("asterviz_trivial." + gr.gid + ".dot", gene_start_end + tp2v_to_string());
	}

	return true;
}
*/

/* remove abutting edge, then call divide_conquer(source, target, res) again */
bool aster::divide_conquer_abutting(int source, int target, aster_result& res)
{	
	assert(source < tp2v.size() && target < tp2v.size());
	int s = tp2v[source];
	int t = tp2v[target];
	assert(s < gr.num_vertices() && t < gr.num_vertices() && s >= 0 && t >= 0);
	assert(s < t - 1);
	assert(source < target - 1);
	assert(gr.out_degree(s) >= 1 || gr.in_degree(t) >= 1);

	if(! gr.edge_exists(s,t)) return false;

	if(verbose >= 2) 
	{
		string msg = "aster processing subgraph, vertex [" + to_string(s) + ", " + to_string(t) + "]"; 
		msg += " (topoIndex [" + to_string(source) + "," + to_string(target) + "]), ";
		msg += "removing the direct edge between them."; 
		cout << msg << endl;
	}

	// remove abutting edge & recusion
	edge_descriptor e = gr.edge(s, t).first;
	assert(e != null_edge);
	double w = gr.get_edge_weight(e);
	assert(gr.edge(e));
	gr.remove_edge(e);
	assert(gr.check_path(s, t));
	divide_conquer(source, target, res);

	int shortestPathIndex = find_shortest_path(res);
	assert(shortestPathIndex >= 0);
	int shortestPathSize = res.subpaths.at(shortestPathIndex).v.size();
	int eventSize        = shortestPathSize - 2;
	assert(eventSize >= 1);
	path p({s, t}, w);
	res.subpaths.push_back(p);
	res.dist += event_size_penalty(eventSize);

	replace_aster_index_to_one_edge(source, target, w);

	if(verbose >= 2) 
	{
		string msg = "aster processed subgraph, vertex [" + to_string(s) + ", " + to_string(t) + "]"; 
		msg += " (topoIndex [" + to_string(source) + "," + to_string(target) + "]), ";
		msg += "with a direct edge between them."; 
		cout << msg << endl;
	}

	return true;
}	

//FIXME: use connected-component to do a more comprehensive work. Now works.
bool aster::divide_conquer_cut_termini(int source, int target, aster_result& res)
{
	assert(source < tp2v.size() && target < tp2v.size());
	int s = tp2v[source];
	int t = tp2v[target];
	assert(s < gr.num_vertices() && t < gr.num_vertices() && s >= 0 && t >= 0);
	assert(s < t - 1);
	assert(source < target - 1);
	assert(gr.out_degree(s) >= 1 || gr.in_degree(t) >= 1);
	if(gr.edge_exists(s,t)) return false;
	if(gr.out_degree(s) == 1 || gr.in_degree(t) == 1) return false;

	// examine if disjoint subgraphs; as the graph is DFS topo-sorted
	vector<pair<int, int>> subgraph_intervals;
	int subgraphNum = divide_conquer_cut_termini_find(source, target, subgraph_intervals);
	if(subgraphNum <= 0) return false;
	if(subgraph_intervals.size() == 0) return false;
	// assertions 
	for(const auto& pp: subgraph_intervals)	
	{
		int subsource = pp.first;
		int subtarget = pp.second;
		int ss = tp2v[subsource];
		int tt = tp2v[subtarget];
		for(int i = subsource; i < subtarget; i++)
		{
			int vidx = tp2v.at(i);
			PEEI peei = gr.out_edges(vidx);
			for(edge_iterator it1 = peei.first, it2 = peei.second; it1 != it2; it1++)
			{
				assert(v2tp.at((*it1)->target()) <= subtarget);
			}
		}
		for(int i = subtarget; i > subsource; i--)
		{
			int vidx = tp2v.at(i);
			PEEI peei = gr.in_edges(vidx);
			for(edge_iterator it1 = peei.first, it2 = peei.second; it1 != it2; it1++)
			{
				assert(v2tp.at((*it1)->source()) >= subsource);
			}
		}
	}

	if(verbose >= 2) 
	{
		string msg = "aster processing subgraph, vertex [" + to_string(s) + ", " + to_string(t) + "]"; 
		msg += " (topoIndex [" + to_string(source) + "," + to_string(target) + "]), ";
		msg += "splitting disjoint graphs at termini"; 
		msg += "\n\tsplitted intervals:";
		for(const auto& iv: subgraph_intervals)
		{
			msg += "\n\t";
			msg += "[" + to_string(tp2v.at(iv.first)) + "-" + to_string(tp2v.at(iv.second)) + "] ";
			msg += "(Topo[" + to_string(iv.first) + "-" + to_string(iv.second) + "]); ";
		}
		cout << msg << endl;
	}

	// get paths for subgraphs
	vector<aster_result> resOfSubgraphs(subgraph_intervals.size());
	for(int i = 0; i < subgraph_intervals.size(); i++)
	{
		aster_result& res1 = resOfSubgraphs[i];
		int subsource = subgraph_intervals[i].first;
		int subtarget = subgraph_intervals[i].second;
		divide_conquer(subsource, subtarget, res1);
		assert(res1.subpaths.size() > 0);
		for(path& p: res1.subpaths)
		{
			assert(p.v.size() >= 1);
			assert(p.v.front() != s);
			assert(p.v.back() != t);
			p.v.insert(p.v.begin(), s);
			p.v.push_back(t);
			assert(origr.valid_path(p.v));
			assert(p.v.size() >= 3);
		}
	}

	// get penalty 
	res.clear();
	int numEvents = 0;
	int sumPenalty = 0;
	for(int i = 0; i < subgraph_intervals.size(); i++)
	{
		aster_result& res1 = resOfSubgraphs[i];
		int shortestPathSize1 = res1.subpaths.at(find_shortest_path(res1)).v.size();
		assert(shortestPathSize1 - 2 >= 1);
		numEvents += shortestPathSize1 - 2;
		sumPenalty += res1.dist;
		res.subpaths.insert(res.subpaths.begin(), res1.subpaths.begin(), res1.subpaths.end());
	}
	assert(numEvents >= 1);
	res.dist = event_size_penalty(numEvents) + sumPenalty;

	double w = 0;
	for(const path& p: res.subpaths) w += p.abd;
	replace_aster_index_to_one_edge(source, target, w);

	if(verbose >= 2) 
	{
		string msg = "aster processed subgraph, vertex [" + to_string(s) + ", " + to_string(t) + "]"; 
		msg += " (topoIndex [" + to_string(source) + "," + to_string(target) + "]), ";
		msg += "with " + to_string(subgraph_intervals.size()) +" disjoint graphs at termini"; 
		cout << msg << endl;
	}

	return true;
}

/* return size of interval pairs 
*  return -1 if not exists 
*  if source/target is cut vertex, 
* 		then all subgraph's sources must be source-vertex's out-edge targets
* 		then all subgraph's targets must be target-vertex's in-edge  sources
*/
int aster::divide_conquer_cut_termini_find(int source, int target, vector<pair<int, int>>& intervals)
{
	assert(source < tp2v.size() && target < tp2v.size());
	int s = tp2v[source];
	int t = tp2v[target];
	assert(s < gr.num_vertices() && t < gr.num_vertices() && s >= 0 && t >= 0);
	assert(s < t - 1);
	assert(source < target - 1);
	assert(gr.out_degree(s) >= 1 || gr.in_degree(t) >= 1);
	if(gr.edge_exists(s,t)) return -1;
	if(gr.out_degree(s) == 1 || gr.in_degree(t) == 1) return -1;

	
	constexpr int IS_SOURCE = 1;
	constexpr int IS_TARGET = 2;
	vector<pair <int, int>> subgraphDisjoinPoints;
	PEEI peei1 = gr.out_edges(s);
	for(edge_iterator it1 = peei1.first, it2 = peei1.second; it1 != it2; it1++)
	{
		edge_descriptor e = *it1;
		int subs = e->target();
		assert(subs != t);
		int subsource = v2tp.at(subs);
		subgraphDisjoinPoints.push_back({subsource, IS_SOURCE});
	}
	PEEI peei2 = gr.in_edges(t);
	for(edge_iterator it1 = peei2.first, it2 = peei2.second; it1 != it2; it1++)
	{
		edge_descriptor e = *it1;
		int subt = e->source();
		assert(subt != s);
		int subtarget = v2tp.at(subt);
		subgraphDisjoinPoints.push_back({subtarget, IS_TARGET});
	}
	sort(subgraphDisjoinPoints.begin(), subgraphDisjoinPoints.end());

	// delete key if 1. IS_SOURCE and previous IS_SOURCE; 2. IS_TARGET and next IS_TARGET
	assert((subgraphDisjoinPoints.begin()->second & IS_SOURCE) >= 1);
	assert((prev(subgraphDisjoinPoints.end())->second & IS_TARGET) >= 1);
    for (auto it = next(subgraphDisjoinPoints.begin()); it != prev(subgraphDisjoinPoints.end()); /*increment separately handled*/ ) 
	{
		if(it->first == next(it)->first) assert((it->second & IS_SOURCE) && (next(it)->second & IS_TARGET));
		int label = it->second;
		assert(label == IS_SOURCE || label == IS_TARGET);
		int prevLable = prev(it)->second;
		int nextLabel = next(it)->second;
        if ((label & IS_SOURCE) && (prevLable & IS_SOURCE)) 
		{
			it->second -= IS_SOURCE;
			if(it->second <= 0) it = subgraphDisjoinPoints.erase(it);
			continue;            
        } 
		else if((label & IS_TARGET) && (nextLabel & IS_TARGET))
		{
			it->second -= IS_TARGET;
			if(it->second <= 0) it = subgraphDisjoinPoints.erase(it);
			continue;
		}
		else 
		{
            ++it;
        }
    }

	// pupulate `intervals`
	// use a vector and sort vector
	intervals.clear();
	for (auto it = subgraphDisjoinPoints.begin(); it != subgraphDisjoinPoints.end(); it++) 
	{
		// multi vertex interval
		int label = it->second;
		assert(label == IS_SOURCE || label == IS_TARGET);
		if(it->second & IS_SOURCE)
		{
			assert(next(it)->second & IS_TARGET);
			int subsource = it->first;
			int subtarget = next(it)->first;
			
			assert(subsource <= subtarget);
			assert(subsource > source);
			assert(subtarget < target);
			intervals.push_back({subsource, subtarget});
		}
		else 
		{
			assert(it->second & IS_TARGET);
			assert(next(it) == subgraphDisjoinPoints.end() || next(it)->second & IS_SOURCE);
		}
    }


	// assertions and validations
	for(const auto& pp: intervals)
	{
		int subsource = pp.first;
		int subtarget = pp.second;
		int ss = tp2v[subsource];
		int tt = tp2v[subtarget];
		for(int i = subsource; i < subtarget; i++)
		{
			int vidx = tp2v.at(i);
			PEEI peei = gr.out_edges(vidx);
			for(edge_iterator it1 = peei.first, it2 = peei.second; it1 != it2; it1++)
			{
				if(v2tp.at((*it1)->target()) <= subtarget) continue;
				intervals.clear();
				return -1;
			}
		}
		for(int i = subtarget; i > subsource; i--)
		{
			int vidx = tp2v.at(i);
			PEEI peei = gr.in_edges(vidx);
			for(edge_iterator it1 = peei.first, it2 = peei.second; it1 != it2; it1++)
			{
				if(v2tp.at((*it1)->source()) >= subsource) continue;
				intervals.clear();
				return -1;
			}
		}
	}

	if(intervals.size() == 0) return -1;
	return intervals.size();
}

/* check if the graph can be divided to two disjoint graphs at a pivot k
** divide and conquer the problem to dnc[s, k], dnc [k, t]
** not disjoint: (s,t) edge exists, or multiple subgraphs can be found
*/
bool aster::divide_conquer_articulation_point(aster_index ai)
{
	int s = ai.s();
	int t = ai.t();
	assert(s < gr.num_vertices() && t < gr.num_vertices() && s >= 0 && t >= 0);
	assert(s < t - 1);
	assert(gr.out_degree(s) >= 1 || gr.in_degree(t) >= 1);
	if(gr.edge_exists(s,t)) return false;
	
	aster_index aileft, airight;
	int pivot = divide_conquer_articulation_find(ai, aileft, airight);
	if (pivot < 0) return false;
	assert(pivot > s && pivot < t);

	divide_conquer(aileft); 	
	PEB peb1 = gr.edge(s, pivot);
	edge_descriptor e1 = peb1.first;
	assert(peb1.second);
	assert(gr.compute_num_paths(s, pivot) == 1);

	divide_conquer(airight);
	PEB peb2 = gr.edge(pivot, t);
	edge_descriptor e2 = peb2.first;
	assert(peb2.second);	
	assert(gr.compute_num_paths(pivot, t) == 1);
	assert(gr.compute_num_paths(s, t) == 1);

	aster_result & res1 = edgeres.at(e1);
	aster_result & res2 = edgeres.at(e2);
	aster_result comb;
	bool combineSuccess = res_combine_consecutive(res1, res2, comb);
	if (!combineSuccess) return false;

	double w = 0;
	for(const path& p: comb.subpaths) w += p.abd;
	replace_aster_index_to_one_edge(ai, w, comb);

	if(verbose >= 2)
	{
		string msg = "aster processed subgraph, vertex [" + to_string(s) + ", " + to_string(t) + "]"; 
		// msg += " (topoIndex [" + to_string(source) + "," + to_string(target) + "]), ";
		msg += "with subgraphs at articulation point " + to_string(tp2v[pivot]);
		msg += " (topoIndex = " + to_string(pivot) + ")";
		cout << msg << endl;
	}

	return true;
}

// find a pivot-TopoIndex s.t. removing this vertex will split grpah to two parts between [source, pivot] and [pivot, target]
// return -1 if cannot fine artivulation point
int aster::divide_conquer_articulation_find(aster_index ai, aster_index left, aster_index right) 
{
	int s = ai.s();
	int t = ai.t();
	if(ai.size() <= 2) return false;
	if(! gr.check_path(s, t)) return false;

	int pivot = -1;
	int index = -1;
	for(int i = 1; i < ai.size() - 1; i++)
	{
		int v = ai.at(i);
		splice_graph gr2(gr);
		gr2.clear_vertex(v);
		if(gr2.check_path(s, t)) continue;
		pivot = v;
		index = i;
		break;
	}
	if (pivot >= t || pivot <= s || pivot <= 0) return -1;

	if(verbose >= 2) 
	{
		string msg = "\t found articulation point = " + to_string(pivot) + ")";
		cout << msg << endl;
	}

	// split ai
	ai.split(index, left, right);
	assert(left.t() == pivot);
	assert(right.s() == pivot);

	//assertion
	splice_graph gr2(gr);
	gr2.clear_vertex(pivot);
	assert(! gr2.check_path(s, t)); 

	return pivot;
}

/*  combines subpaths of left and right sides of articulation point k
*	res1, res2 could be empty
*/
bool aster::res_combine_consecutive(aster_result& res1,  aster_result& res2, aster_result& comb) const
{
	if(res1.subpaths.size() == 0 && res2.subpaths.size() == 0) return false;
	if(res1.subpaths.size() == 0) 
	{
		comb = res2;
		return true;
	}
	if(res2.subpaths.size() == 0) 
	{
		comb = res1;
		return true;
	}

	assert(res1.subpaths.size() > 0);
	assert(res1.subpaths[0].v.size() > 0);
	assert(res2.subpaths.size() > 0);
	assert(res2.subpaths[0].v.size() > 0);

	int k = res1.subpaths[0].v.back();

	for(const path& p: res1.subpaths)	
	{
		assert(p.v.size() >= 1); 
		assert(p.v.back() == k);
	}
	for(const path& p: res2.subpaths)	
	{
		assert(p.v.size() >= 1); 
		assert(p.v.front() == k);
	}

	int index1 = find_shortest_path(res1);
	int index2 = find_shortest_path(res2);
	assert(index1 >= 0);
	assert(index2 >= 0);
	assert(index1 < res1.subpaths.size());
	assert(index2 < res2.subpaths.size());
	
	comb.clear();
	const path& rAnchor = res2.subpaths[index2];
	for(int i = 0; i < res1.subpaths.size(); i++)	
	{
		if(i == index1) continue;
		const path& p = res1.subpaths.at(i);
		double abd = p.abd < rAnchor.abd? p.abd: rAnchor.abd;
		vector<int> v;
		v.insert(v.end(), p.v.begin(), p.v.end());
		v.insert(v.end(), next(rAnchor.v.begin()), rAnchor.v.end());
		assert(origr.valid_path(v));
		comb.subpaths.push_back(path(v, abd));
	}
	const path& lAnchor = res1.subpaths.at(index1);
	for(int i = 0; i < res2.subpaths.size(); i++)	
	{
		if(i == index2) continue;
		const path& p = res2.subpaths.at(i);
		double abd = p.abd < lAnchor.abd? p.abd: lAnchor.abd;
		vector<int> v;
		v.insert(v.end(), lAnchor.v.begin(), lAnchor.v.end());
		v.insert(v.end(), next(p.v.begin()), p.v.end());
		assert(origr.valid_path(v));
		comb.subpaths.push_back(path(v, abd));
	}
	double abd = lAnchor.abd < rAnchor.abd? lAnchor.abd: rAnchor.abd;
	vector<int> v;
	v.insert(v.end(), lAnchor.v.begin(), lAnchor.v.end());
	v.insert(v.end(), next(rAnchor.v.begin()), rAnchor.v.end());
	assert(origr.valid_path(v));
	comb.subpaths.push_back(path(v, abd));

	assert(comb.subpaths.size() >= res1.subpaths.size() + res2.subpaths.size() - 1);

	if(verbose >= 3)
	{
		cout << "combining subpaths paths squeezing articulation point " << k << endl ;
		cout << "paths on the left side of art-point: ";
		for(const path& p: res1.subpaths)
		{
			cout << endl;
			cout << "weight = " << p.abd << " ; path = [";
			printv(p.v);
			cout << "]";
		}
		cout << endl;
		
		cout << "paths on the right side of art-point: ";
		for(const path& p: res2.subpaths)
		{
			cout << endl;
			cout << "weight = " << p.abd << " ; path = [";
			printv(p.v);
			cout << "]";
		}
		cout << endl;
		
		cout << "combined paths ";
		for(const path& p: comb.subpaths)
		{
			cout << endl;
			cout << "weight = " << p.abd << " ; path = [";
			printv(p.v);
			cout << "]";
		}
		cout << endl;
	}

	return true;
}

bool aster::res_combine_parallel(aster_result& res1,  aster_result& res2, aster_result& comb) const
{	
	if(res1.subpaths.size() == 0 && res2.subpaths.size() == 0) return false;
	if(res1.subpaths.size() == 0) 
	{
		comb = res2;
		return true;
	}
	if(res2.subpaths.size() == 0) 
	{
		comb = res1;
		return true;
	}
	assert(res1.subpaths.size() > 0);
	assert(res1.subpaths[0].v.size() > 0);
	assert(res2.subpaths.size() > 0);
	assert(res2.subpaths[0].v.size() > 0);

	int s = res1.subpaths[0].v.front();
	int t = res1.subpaths[0].v.back();

	for(const path& p: res1.subpaths)	
	{
		assert(p.v.size() >= 1); 
		assert(p.v.front() == s);
		assert(p.v.back() == t);
	}
	for(const path& p: res2.subpaths)	
	{
		assert(p.v.size() >= 1); 
		assert(p.v.front() == s);
		assert(p.v.back() == t);
	}

	comb.clear();
	comb.subpaths.insert(comb.subpaths.end(), res1.subpaths.begin(), res1.subpaths.end());
	comb.subpaths.insert(comb.subpaths.end(), res2.subpaths.begin(), res2.subpaths.end());
	assert(comb.subpaths.size() == res1.subpaths.size() + res2.subpaths.size());

	//TODO: calculate dist
	return true;
}

/* examine if dnc unitig between source to target; if true, populate res */
bool aster::divide_conquer_unitig(aster_index ai)
{
	int s = ai.s();
	int t = ai.t();
	assert(s < gr.num_vertices() && t < gr.num_vertices() && s >= 0 && t >= 0);
	assert(s < t);
	assert(gr.out_degree(s) >= 1 || gr.in_degree(t) >= 1);

	if(verbose >= 3)  //CLEAN:
	{
		string msg = "aster checking unitig in subgraph, vertex [" + to_string(s) + ", " + to_string(t) + "]"; 
		cout << msg;
		cout << endl;
	}

	aster_result unitigRes;
	bool   _avg_ = false;       							// average if true, geom mean if false
	double w     = _avg_? 0: 1;

	int n = gr.compute_num_paths(s, t, 2);
	if(s == t - 1)  assert(n == 1); // only two vertices
	assert(n >= 1);
	if(n > 1) return false;	

	// populate unitig
	double c = 0.0;
	int ss = s;
	while(true)
	{
		assert(ss <= t);
		if (ss == t) break;
		assert(gr.out_degree(ss) >= 1);
		if (gr.out_degree(ss) > 1) return false;
		edge_descriptor e = (*(gr.out_edges(ss).first));
		double ew = gr.get_edge_weight(e);
		w = _avg_? (w + ew): (w * ew);
		ss  = e->target();

		aster_result __res__;
		res_combine_consecutive(unitigRes, edgeres.at(e), __res__);
		unitigRes = __res__;
		c += 1.0;
	}
	assert(valid_paths(unitigRes));
	w = _avg_? (w / c): pow(w, 1.0/ c);

	aster_result __res__;
	res_combine_consecutive(unitigRes, edgeres.at(e), __res__);
	unitigRes = std::move(__res__);		//FIXME:TODO: is this correct?	

	if(verbose >= 2) 
	{
		string msg = "aster processed subgraph, vertex [" + to_string(s) + ", " + to_string(t) + "]"; 
		// msg += " (topoIndex [" + to_string(source) + "," + to_string(target) + "]), ";
		msg += "with a unitig path: "; 
		cout << msg;
		cout << endl;
	}

	replace_aster_index_to_one_edge(ai, w, unitigRes);
	return true;
}

/* examine if dnc single vertex; if true, populate res 
*/
bool aster::divide_conquer_single_vertex(aster_index ai)
{
	int s = ai.s();
	int t = ai.t();
	assert(s < gr.num_vertices() && t < gr.num_vertices() && s >= 0 && t >= 0);
	assert(s <= t);
	if(s != t) return false;
	
	if(verbose >= 2) 
	{
		string msg = "aster processed subgraph, vertex [" + to_string(s) + ", " + to_string(t) + "]"; 
		// msg += " (topoIndex [" + to_string(source) + "," + to_string(target) + "]), ";
		msg += "with a single-vertex [" + to_string(s) + ", " + to_string(t) + "]"; 
		cout << msg << endl;
	}

	if(gr.degree(s) == 0)  return true;

	aster_result res;
	double abd = gr.get_vertex_weight(s);

	throw runtime_error("single vertex should not be called");
	return false;
}

int aster::init_edgeres()
{
	edgeres.clear();
	PEEI pei = gr.edges();
	for(edge_iterator it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
	{
		edge_descriptor e = *it1;
		double w = gr.get_edge_weight(e);
		edgeres.insert({e, aster_result(vector<int>{e->source(), e->target()}, w)});
	}
	return 0;
}

/*
* DFS based topological sort 
* This guarantees all disjoint subgraphs are gathered together
* assuming no two vertices have the same <lpos, rpos>
*/
int aster::topological_sort_vertices()
{
	assert(gr.num_vertices() >= 2);
	tp2v.clear();
	v2tp.clear();

	vector<bool> visited(gr.num_vertices(), false);
	for(int i = 0; i < gr.num_vertices(); i++)	
	{
		topological_sort_vertices_visit(i, visited, tp2v);
	}
	reverse(tp2v.begin(), tp2v.end());
	v2tp.resize(tp2v.size());
	for(int i = 0; i < tp2v.size(); i++) v2tp[tp2v[i]] = i;

	/* // assertions topo sorted
	for(int i = 0; i < gr.num_vertices(); i++)	
	{
		if(i >= 1) assert(gr.get_vertex_info(i - 1).rpos <= gr.get_vertex_info(i).lpos);
		if (i < gr.num_vertices() - 1)  // proceeding out edges
		{
			PEEI          ei  = gr.out_edges(i);
			edge_iterator it1 = ei.first;
			edge_iterator it2 = ei.second;
			for(; it1 != it2; it1++) assert(gr.get_vertex_info(i).rpos <= gr.get_vertex_info((*it1)->target()).lpos);
		}
		if (i >= 1)  // trailing in edges
		{
			PEEI          ei  = gr.in_edges(i);
			edge_iterator it1 = ei.first;
			edge_iterator it2 = ei.second;
			for(; it1 != it2; it1++) assert(gr.get_vertex_info(i).lpos >= gr.get_vertex_info((*it1)->source()).rpos);
		}
	} */

	if(verbose >= 3)	cout << "aster sorted " << tp2v_to_string() << endl;
	assert(tp2v.front() == 0);
	assert(tp2v.back() == gr.num_vertices() - 1);
	assert(tp2v.size() == gr.num_vertices());
	
	return 0;
}

/* visit a node i in DFS, push i to _tp2v */
int aster::topological_sort_vertices_visit(int i, vector<bool>& visited, vector<int>& _tp2v)
{
	if(visited[i]) return 0;
	visited[i] = true;

	if(gr.out_degree(i) == 0)
	{
		_tp2v.push_back(i);
		return 0;
	}

	PEEI peei;
	edge_iterator it1, it2;
	for(peei = gr.out_edges(i), it1 = peei.first, it2 = peei.second; it1 != it2; it1++)
	{
		int j = (*it1)->target();
		topological_sort_vertices_visit(j , visited, _tp2v);
	}
	_tp2v.push_back(i);
	return 0;
}

/* Edges are sorted by their <source, target> represented by tp2v
*  edges should be sorted & fetched independent of splice_graph implementation
*/ 
int aster::topological_sort_index_edges()
{
	MEI e2i;							// edge map, from edge to index, sorted by position
	VE i2e;								// edge map, from index to edge, sorted	by position
	map<pair<int, int>, edge_descriptor> sortedEdges;
	PEEI pei = gr.edges(); 
	for(edge_iterator it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
	{
		edge_descriptor e = *it1;
		sortedEdges.insert({{v2tp[e->source()], v2tp[e->target()]}, e});	
	}
	for(auto it = next(sortedEdges.begin()); it != sortedEdges.end(); it++) 
	{
		auto pv = prev(it); 
		int p1 = it->first.first, p2 = it->first.second;
		int q1 = pv->first.first, q2 = pv->first.second;
		assert(p1 >= q1);
		assert(p1 > q1 || (p1 == q1 && p2 > q2));
	}

	i2e.clear();
	e2i.clear();
	int index = 0;
	for(const auto& iie: sortedEdges)
	{
		edge_descriptor e = iie.second;
		e2i.insert({e, index});
		i2e.push_back(e);
		index++;
	}
	assert(e2i.size() == index); 
	assert(i2e.size() == index);
	return 0;
}

/*
*  aggresively remove intersecting edges, whichever is topilocially smaller
*/
/* bool aster::aggressive_purge_intersecting_edges()
{
	if(gr.num_vertices() > 1000)
	{
		cout << "graph too big to pruge, #vertex = " << gr.num_vertices() << endl;
		return 0;
	}

	int purgeCount = 0;
	for (int i = 0; i < i2e.size() - 1; i ++)
	{
		edge_descriptor edge1 = i2e[i];
		if (edge1 == null_edge) continue;
		for (int j = i + 1; j < i2e.size(); j ++)
		{
			edge_descriptor edge2 = i2e[j];
			if (edge2 == null_edge) continue;
			if (!gr.intersect(edge1, edge2)) continue;
			i2e[i] = null_edge;
			e2i.erase(edge1);
			gr.remove_edge(edge1);
			purgeCount ++;
			break;
		}
	}
	// if (verbose >= 2 && purgeCount <  1) cout << "graph does not have intersecting edges." << endl;
	if (verbose >= 2 && purgeCount >= 1) cout << "removed " << purgeCount << " intersecting edges from graph.\nSort graph again" << endl;
	gr.refine_splice_graph();

	if (purgeCount >= 1)
	{
		topological_sort_vertices();
		topological_sort_index_edges();
	}
	return (purgeCount >= 1);
} */

int aster::balance_vertex(int vertexIndex)
{
	int v = vertexIndex;
	if(gr.degree(v) <= 0) return 0;

	edge_iterator it1, it2;
	PEEI pei;
	double w1 = 0, w2 = 0;
	for(pei = gr.in_edges(v), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
	{
		double w = gr.get_edge_weight(*it1);
		w1 += w;
	}
	for(pei = gr.out_edges(v), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
	{
		double w = gr.get_edge_weight(*it1);
		w2 += w;
	}

	assert(w1 >= SMIN);
	assert(w2 >= SMIN);

	// use sqrt-meature
	double ww = sqrt(w1 * w2);

	double r1 = ww / w1;
	double r2 = ww / w2;

	double m1 = 0, m2 = 0;
	for(pei = gr.in_edges(v), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
	{
		double wx = gr.get_edge_weight(*it1);
		double wy = wx * r1;
		if(wy < min_guaranteed_edge_weight)
		{
			m1 += (min_guaranteed_edge_weight - wy);
			wy = min_guaranteed_edge_weight;
		}
		gr.set_edge_weight(*it1, wy);
	}
	for(pei = gr.out_edges(v), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
	{
		double wx = gr.get_edge_weight(*it1);
		double wy = wx * r2;
		if(wy < min_guaranteed_edge_weight)
		{
			m2 += min_guaranteed_edge_weight - wy;
			wy = min_guaranteed_edge_weight;
		}
		gr.set_edge_weight(*it1, wy);
	}

	if(m1 > m2)
	{
		edge_descriptor e = gr.max_out_edge(v);
		double w = gr.get_edge_weight(e);
		gr.set_edge_weight(e, w + m1 - m2);
	}
	else if(m1 < m2)
	{
		edge_descriptor e = gr.max_in_edge(v);
		double w = gr.get_edge_weight(e);
		gr.set_edge_weight(e, w + m2 - m1);
	}

	return 0;
}

int aster::remove_small_junctions()
{
	SE se;
	for(int i = 1; i < gr.num_vertices() - 1; i++)
	{
		if(gr.degree(i) <= 0) continue;

		bool b = true;
		edge_iterator it1, it2;
		PEEI pei;
		int32_t p1 = gr.get_vertex_info(i).lpos;
		int32_t p2 = gr.get_vertex_info(i).rpos;
		double wi = gr.get_vertex_weight(i);

		// compute max in-adjacent edge
		double ws = 0;
		for(pei = gr.in_edges(i), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
		{
			edge_descriptor e = (*it1);
			int s = e->source();
			double w = gr.get_vertex_weight(s);
			if(s == 0) continue;
			if(gr.get_vertex_info(s).rpos != p1) continue;
			if(w < ws) continue;
			ws = w;
		}

		// remove small in-junction
		for(pei = gr.in_edges(i), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
		{
			edge_descriptor e = (*it1);
			int s = e->source();
			double w = gr.get_edge_weight(e);
			if(s == 0) continue;
			if(gr.get_vertex_info(s).rpos == p1) continue;
			if(ws < 2.0 * w * w + 18.0) continue;
			if(wi < 2.0 * w * w + 18.0) continue;

			se.insert(e);
		}

		// compute max out-adjacent edge
		double wt = 0;
		for(pei = gr.out_edges(i), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
		{
			edge_descriptor e = (*it1);
			int t = e->target();
			double w = gr.get_vertex_weight(t);
			if(t == gr.num_vertices() - 1) continue;
			if(gr.get_vertex_info(t).lpos != p2) continue;
			if(w < wt) continue;
			wt = w;
		}

		// remove small in-junction
		for(pei = gr.out_edges(i), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
		{
			edge_descriptor e = (*it1);
			double w = gr.get_edge_weight(e);
			int t = e->target();
			if(t == gr.num_vertices() - 1) continue;
			if(gr.get_vertex_info(t).lpos == p2) continue;
			if(wt < 2.0 * w * w + 18.0) continue;
			if(wi < 2.0 * w * w + 18.0) continue;

			se.insert(e);
		}

	}

	if(se.size() <= 0) return false;

	for(SE::iterator it = se.begin(); it != se.end(); it++)
	{
		edge_descriptor e = (*it);
		vertex_info v1 = gr.get_vertex_info(e->source());
		vertex_info v2 = gr.get_vertex_info(e->target());
		if(verbose >= 2) printf("remove small junction: length = %d, pos = %d-%d\n", v2.lpos - v1.rpos, v2.lpos, v1.rpos);
		gr.remove_edge(e);
	}

	return true;
}

int aster::edge_path_to_vertex_path(const VE& edgePath, VI& vertexPath) const
{
	if (edgePath.size() == 0) return 0;
	vertexPath.clear();
	vertexPath.resize(edgePath.size() + 1);
	vertexPath[0] = edgePath.front()->source();
	int i = 1;
	for(const edge_descriptor& e: edgePath)
	{
		int s = e->source();
		int t = e->target();
		assert(s == vertexPath[i - 1]);
		vertexPath[i] = t;
		i ++;
	}
	assert(vertexPath.size() >= 2);
	assert(origr.valid_path(vertexPath));
	return 0;
}

// return index of longest path in res.subpaths
int aster::find_longest_path(const aster_result& res) const
{
	int longestPathIndex = -1;
	int longestPathSize = -1;
	if (res.subpaths.size() == 0) return longestPathIndex;

	longestPathIndex = 0;
	longestPathSize = res.subpaths.front().v.size();
	for(int i = 0; i < res.subpaths.size(); i++) 
	{
		const path& p = res.subpaths.at(i);
		if(p.v.size() <= longestPathSize) continue;
		longestPathSize = p.v.size();
		longestPathIndex = i;
	}
	return longestPathIndex;
}

// return index of shortest path in res.subpaths
int aster::find_shortest_path(const aster_result& res) const
{
	int shortestPathIndex = -1;
	int shortestPathSize = -1;
	if (res.subpaths.size() == 0) return shortestPathIndex;

	shortestPathIndex = 0;
	shortestPathSize = res.subpaths.front().v.size();
	for(int i = 0; i < res.subpaths.size(); i++) 
	{
		const path& p = res.subpaths[i];
		if(p.v.size() >= shortestPathSize) continue;
		shortestPathSize = p.v.size();
		shortestPathIndex = i;
	}
	return shortestPathIndex;
}

/* 
*  Assume internal nodes must be "closed" in [s, t]
*  only remove edges in [s,t] (toposorted). 
*  Maybe not necessarily removes all edges from/to vertices in [s,t] if it is not closed
*  i.e. for nodes reachable from s and reachable to t, they cannot reach any other nodes but themselves and s,t
*/
edge_descriptor aster::replace_aster_index_to_one_edge(aster_index ai, double w, aster_result res)
{
	int s = ai.s();
	int t = ai.t();
	int source = v2tp.at(s);
	int target = v2tp.at(t);
	assert(s < gr.num_vertices() && t < gr.num_vertices() && s >= 0 && t >= 0);
	assert(s < t);
	assert(gr.out_degree(s) >= 1 || gr.in_degree(t) >= 1);	

	// under 'closed' assumption all edges in (s, t) open interval must have sources/targets in [s, t] closed interval
	for(int k: ai.get_index())
	{
		if(gr.degree(k) <= 0) continue;
		assert(k >= s);
		assert(k <= t);
		set<edge_descriptor> setEdge;
		PEEI pi = gr.in_edges(k);
		for(edge_iterator it = pi.first; it != pi.second; it++)
		{
			setEdge.insert(*it);
		}
		PEEI po = gr.out_edges(k);
		for(edge_iterator it = po.first; it != po.second; it++)
		{
			setEdge.insert(*it);
		}
		for(edge_descriptor e: setEdge)
		{
			assert(gr.edge(e));
			if(! ai.find_index(e->source())) continue;;
			if(! ai.find_index(e->target())) continue;
			gr.remove_edge(e);
			edgeres.erase(e);
		}
	}
	assert(! gr.check_path(s, t));

	// put edge & res
	edge_descriptor e_new = gr.add_edge(s, t);
	edge_info ei;
	ei.weight = w;
	gr.set_edge_info(e_new, ei);
	gr.set_edge_weight(e_new, w);
	assert(e_new != NULL);
	assert(gr.ewrt.find(e_new) != gr.ewrt.end());
	assert(gr.einf.find(e_new) != gr.einf.end());
	edgeres.insert({e_new, res});

	assert(!gr.refine_splice_graph());
	return e_new;
}

bool aster::valid_paths(aster_result res) const
{
	return valid_paths(res.subpaths);
}

bool aster::valid_paths(vector<path> paths) const
{
	for(const path& p : paths)
	{
		if (!origr.valid_path(p.v)) return false;
	}
	return true;
}

// assign path.nf, populate trsts and non_full_trsts
int aster::get_transcripts()
{
	if(successStatus == false) return 0;
	if(mode == aster_mode::STAT_ONLY) return 0;
	if(origr.num_edges() == 0) return 0;
	if(origr.num_vertices() == 2) return 0;
	if(paths.size() == 0) return 0;

	//assign nf
	for(path& p : paths)
	{
		vector<int>& v = p.v;
		bool empty = false;
		for(int i = 0; i < v.size(); i++)
		{
			if(gr.get_vertex_info(v[i]).type == EMPTY_VERTEX) 
			{
				empty = true;
				break;
			}
		}
		p.nf = empty? 1:0;
	}

	//validate s-t path
	for(const path& p : paths)
	{
		const vector<int>& v = p.v;
		assert(v.front() == 0);
		assert(v.back() == origr.num_vertices() - 1);
		assert(origr.valid_path(v));
	}

	trsts.clear();	
	non_full_trsts.clear();
	origr.output_transcripts1(trsts, non_full_trsts, paths);
	return 0;
}


int aster::make_stats()
{
	if(verbose >= 3 && num_graph != 0 && num_graph % 100 == 0)	print_stats();

	num_graph ++;
	for(int i = 1; i < gr.num_vertices() - 1; i++)
	{
		if(gr.degree(i) == 0) continue;
		num_exon ++;
	}

	PEEI pei = gr.edges();
	for (auto it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
	{
		edge_descriptor e = *it1;
		if(e->source() == 0) continue; 
		if(e->target() == gr.num_vertices() - 1) continue;
		num_intron++;
	}	

	return 0;
}

int aster::print_stats()
{
	cout << "aster stats ================================================================" << endl;
	cout << "\t num graphs = " << num_graph  << endl; 
	cout << "\t num intersecting graphs " << num_intersecting_graph << endl;
	cout << "\t num exons = " << num_exon  << endl; 
	cout << "\t num introns = " << num_intron << endl;
	cout << "\t num intersecting introns = " << num_intersecting_intron_count << endl;
	cout << "\t num intersecting intron pairs " << num_intersecting_intron_pair << endl;
	cout << "\t D&C algorithm counter: ";
	cout << "\t " << dnc_counter_single;
	cout << "\t " << dnc_counter_unitig;
	cout << "\t " << dnc_counter_abutting;
	cout << "\t " << dnc_counter_cut_vertex;
	cout << "\t " << dnc_counter_articulation_point_disjoint;
	cout << "\t " << dnc_counter_resolve_trivial_intersection;
	cout << endl;
	cout << "============================================================================" << endl;
	return 0;
}

string aster::tp2v_to_string() const
{
	string tp2vString = gr.gid + "\n\tDFS TopoSorted vertex index vector:";
	for (int i = 0; i < tp2v.size(); i++)
	{
		if (i % 10 == 0)  tp2vString += "\n\t[" + to_string(i) + "]:";
		tp2vString += to_string(tp2v[i]);
		tp2vString += " ";
	}
	return tp2vString;
}

// exponential penalty guarantees to violate triangle inequality
int aster::event_size_penalty(int eventSize) const
{
	assert(eventSize >= 0);
	return pow(2, eventSize) - 1;
}

/*
* return: 
*		positive int: number of edits
*		assertion error: size not positive
*/
int aster::path_distance(const path& p1, const path& p2) const
{
	throw runtime_error("aster::path_distancen not implemented yet");
	const vector<int>& v1 = p1.v; 
	const vector<int>& v2 = p2.v;
	assert(v1.size() > 0);
	assert(v2.size() > 0);
	if(v1.back() > v2.front()) return -1;
	if(v2.back() > v1.front()) return -1;

	int edits = - basic_algo::ref_sw_query_nw(v1, v2, -1, -2, 0);

	return edits; //TODO: not finished
}

