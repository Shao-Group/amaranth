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

aster::aster(const splice_graph &g, const hyper_set &h, bool _avg)
	: origr(g), gr(g), hs(h), avgMode(_avg)
{
	paths.clear();				// predicted paths, original v index, inclusive
    trsts.clear();				// predicted transcripts, original v index
	non_full_trsts.clear();		// predicted non full length transcripts
	mode = asterMode;
	strategy = asterStrategy;

    topological_sort_vertices();
	prepare_graph();
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

	if(gr.num_vertices() > 1000) //TODO:CLEAN:
	{
		cerr << gr.gid <<" graph is very big to process with D&C, but still do " << gr.num_vertices() << endl; 
		// successStatus = false;
		// return 0;
	}	

	try
	{
		divide_conquer();
		successStatus = true;
	}
	catch(const aster_error& e)
	{
		assert(0); // should not have err
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
	aster_index ai({tp2v});
	divide_conquer(ai);
	
	if(false && output_graphviz_files && verbose >= 3) 
	{
		string gene_start_end = gr.chrm + ":"
							  + to_string(gr.get_vertex_info(0).lpos) + "-"
							  + to_string(gr.get_vertex_info(gr.num_vertices() - 1).rpos)
							  + "\n";
		gr.graphviz("asterviz.final" + gr.gid + ".dot", gene_start_end + tp2v_to_string());
	}
	// collect paths
	int s = 0;
	int t = gr.num_vertices() - 1;
	assert(gr.num_edges() == 1);
	assert(gr.edge_exists(s,t));
	edge_descriptor e = gr.edge(s, t).first;
	paths.clear();
	paths = edgeres.at(e).subpaths;
	if(paths.size() >= 1) for(const path & p : paths)	assert(origr.valid_path(p.v));
	return 0;
}

// make a copy of graph with only edges inside aster_index, to local splice graph
// nodes are not removed
int aster::local_graph(const aster_index& ai, splice_graph& local)
{
	assert_closed_vertex_interval(ai);
	MEE x2y, y2x;
	local.clear();
	local.copy(gr, x2y, y2x);
	for(int i = 0; i < local.num_vertices(); i++) // TODO: optimize
	{
		if(ai.index_found(i)) continue;
		local.clear_vertex(i);
	}
	return 0;
}

// divide_conquer(i ,j) solves a subproblem between 
// i, j are tp2v indices
int aster::divide_conquer(aster_index ai)
{
	stepCount ++;
	int stepCountLocal = stepCount;
	non_isolated_vertex_index(ai);
	assert_closed_vertex_interval(ai);

	int s = ai.s();
	int t = ai.t();
	splice_graph local;
	local_graph(ai, local);
	assert(local.num_vertices() == gr.num_vertices());

	if(verbose >= 2)
	{
		string msg;
		msg += "aster processing subgraph, vertex [" + to_string(s) + ", " + to_string(t) + "] ";
		msg += "(step " + to_string(stepCount) + ")";
		cout << msg << endl;
	}

	if(false && verbose >= 2 && output_graphviz_files)
	{
		string gene_start_end = gr.chrm + ":"
							+ to_string(gr.get_vertex_info(0).lpos) + "-"
							+ to_string(gr.get_vertex_info(gr.num_vertices() - 1).rpos)
							+ "\n";
		gr.graphviz("asterviz." + gr.gid + ".step" + to_string(stepCount) + ".dot", gene_start_end + tp2v_to_string());
	}

	assert(s <= t);
	assert(s < gr.num_vertices() && t < gr.num_vertices() && s >= 0 && t >= 0);

	if (doesFastDnC && resolve_trivial_node(ai))		
	{
		dnc_counter_resolve_trivial_node ++;
		if(verbose >= 3) cout << "aster completed step " << stepCountLocal << endl;
		return 0;
	}
	if (divide_conquer_single_vertex(ai))			
	{
		dnc_counter_single ++;
		if(verbose >= 3) cout << "aster completed step " << stepCountLocal << endl;
		return 0;
	} 
	if (divide_conquer_unitig(ai, local))		
	{
		dnc_counter_unitig ++;
		if(verbose >= 3) cout << "aster completed step " << stepCountLocal << endl;
		return 0;
	}
	if (divide_conquer_abutting(ai))				
	{
		dnc_counter_abutting ++;
		if(verbose >= 3) cout << "aster completed step " << stepCountLocal << endl;
		return 0;
	}
	if (divide_conquer_cut_termini(ai)) 	
	{
		dnc_counter_cut_vertex ++;
		if(verbose >= 3) cout << "aster completed step " << stepCountLocal << endl;
		return 0;
	}
	if (divide_conquer_articulation_point(ai))		
	{
		dnc_counter_articulation_point_disjoint ++;
		if(verbose >= 3) cout << "aster completed step " << stepCountLocal << endl;
		return 0;
	}
	if (resolve_trivial_node(ai))		
	{
		dnc_counter_resolve_trivial_node ++;
		if(verbose >= 3) cout << "aster completed step " << stepCountLocal << endl;
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
		if(verbose >= 3) cout << "aster completed step " << stepCountLocal << endl;
		return 0;
	}
	if(mode != aster_mode::MINI && greedy(ai))
	{
		if(verbose >= 3) cout << "aster completed step " << stepCountLocal << endl;
		return 0;
	}

	num_intersecting_graph ++;

	string msg = "aster-mini failed on graph " + gr.gid;
	msg += " [" + gr.chrm + ":" + to_string(gr.get_vertex_info(0).lpos) + "-" + to_string(gr.get_vertex_info(0).rpos) + "]";
	msg += " in subgraph vertexIndex [" + to_string(s) + ", " + to_string(t) + "]";
	msg += " error at step " + to_string(stepCountLocal);
	
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

/*	find one trivial node and decompose it
	Decomposed node is poped into the edgeres of its decomposed edges
	s and t of the trivial node is connected, original edges removed
*/
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
		if(gr.degree(v) == 0) continue;
		if(gr.in_degree(v) > 1 && gr.out_degree(v) > 1) continue;

		if(verbose >= 3) 
		{
			string msg = "aster processing trivial node " + to_string(v); 
			msg += " in subgraph, vertex [" + to_string(s) + ", " + to_string(t) + "]"; 
			cout << msg << endl;
		}

		PEEI inEdges = gr.in_edges(v);
		PEEI outEdges = gr.out_edges(v);
		assert((gr.in_degree(v) == 1 || gr.out_degree(v) == 1) && gr.degree(v) >= 2);
		edges_combine_consecutive_and_replace(inEdges, outEdges);
		assert(gr.degree(v) == 0);	

		if(verbose >= 2) 
		{
			string msg = "aster resolved a trivial node " + to_string(v); 
			msg += " in subgraph, vertex [" + to_string(s) + ", " + to_string(t) + "]"; 
			cout << msg << endl;
		}
		ai.erase_index(i);
		divide_conquer(ai);
		return true;
	}
	return false;
}

/* 	consecutively combine edges, and replace original edges with newly combined ones
	1. combine their res and push_back to edgeres
	2. remove those edges' edgeres
	2. remove those edges
	Does: combine inEdges res in parallel, combine out Edges res in parallel, then combine In/Out consecutive
*/
int aster::edges_combine_consecutive_and_replace(PEEI inEdges, PEEI outEdges)
{
	//combine res
	aster_result resIn = edgeres.at(*(inEdges.first)) ;
	int inResCount = 0;
	set<edge_descriptor> originalInEdges;
	for(edge_iterator it1 = inEdges.first, it2 = inEdges.second; it1 != it2; it1++)
	{
		edge_descriptor in = *it1;
		originalInEdges.insert(in);
		aster_result resCombIn;
		auto& res1 = edgeres.at(in);
		inResCount += res1.subpaths.size();
		if (it1 != inEdges.first) 
		{
			res_combine_parallel(resIn, res1, resCombIn, false, true);
			resIn = resCombIn;
		}
		assert(resIn.size() >= 1);
	}
	aster_result resOut = edgeres.at(*(outEdges.first));
	int outResCount = 0;
	set<edge_descriptor> originalOutEdges;
	for(edge_iterator ot1 = outEdges.first, ot2 = outEdges.second; ot1 != ot2; ot1++)
	{
		edge_descriptor out = *ot1;
		originalOutEdges.insert(out);	
		aster_result resCombOut;
		auto& res2 = edgeres.at(out);
		outResCount += res2.subpaths.size();
		if(ot1 != outEdges.first) 
		{
			res_combine_parallel(resOut, res2, resCombOut, true, false);
			resOut = resCombOut;
		}
		assert(resOut.size() >= 1);
	}
	aster_result resComb;
	assert(resIn.size() >= 1);
	assert(resOut.size() >= 1);
	res_combine_consecutive(resIn, resOut, resComb);
	assert(resComb.size() >= 1);
	assert(resComb.subpaths.size() == inResCount + outResCount - 1);
	assert(originalInEdges.size() == 1 || originalOutEdges.size() == 1);

	// pop new edgeres from resComb
	map<pair<int, int>, vector<const path*> > st2Path;
	for(const path& p: resComb.subpaths)
	{
		int s = p.v.front();
		int t = p.v.back();
		if(st2Path.find({s,t}) == st2Path.end()) st2Path[{s,t}] = {};
		st2Path[{s,t}].push_back(&p);
	}
	for(edge_descriptor in: originalInEdges)
	{
		for(edge_descriptor out: originalOutEdges)
		{
			edge_combine_consecutive_pop_res(in, out, st2Path);
		}
	}

	// remove
	for(edge_descriptor e: originalInEdges)
	{
		edgeres.erase(e);
		gr.remove_edge(e);
	}
	for(edge_descriptor e: originalOutEdges)
	{
		edgeres.erase(e);
		gr.remove_edge(e);
	}
	return 0;
}


/* 	consecutively combine 2 edges, WITHOUT replacement
	i.e.	combine their res and push_back to edgeres
*/
//FIXME:TODO: hyperset?
int aster::edge_combine_consecutive_pop_res(edge_descriptor in, edge_descriptor out, const map<pair<int, int>, vector<const path*> >& st2Path)
{
	assert(in->target() == out->source());
	int source = in->source();
	int target = out->target();

	// new weight
	double w1 = gr.get_edge_weight(in);
	double w2 = gr.get_edge_weight(out);
	bool   _avg_ = avgMode;	
	double w = _avg_? (w1 + w2 / 2.0) : pow(w1 * w2, 1.0/2.0);

	// new res
	aster_result combinedRes;
	for(const path* pptr :st2Path.at({source, target}))
	{
		combinedRes.subpaths.push_back(*pptr);
	}
	assert(combinedRes.subpaths.size() >= 1);

	// put edge & res
	if(PEB peb = gr.edge(source, target); peb.second)
	{
		edge_descriptor e = peb.first;
		edge_info ei =  gr.get_edge_info(e) ;
		double weight = ei.weight + w;
		ei.weight = weight;
		gr.set_edge_info(e, ei);
		gr.set_edge_weight(e, weight);
		assert(gr.ewrt.find(e) != gr.ewrt.end());
		assert(gr.einf.find(e) != gr.einf.end());
		aster_result res;
		res_combine_parallel(combinedRes, edgeres.at(e), res);
		edgeres[e] = res;
		assert(res.subpaths.size() >= 1);
		for(const path& p: res.subpaths)
		{
			assert(p.v.front() == source);
			assert(p.v.back() == target);
		}
	}
	else
	{
		assert(! gr.edge_exists(source, target));
		edge_descriptor e_new = gr.add_edge(source, target);
		assert(edgeres.find(e_new) == edgeres.end());
		edgeres[e_new] = combinedRes;
		edge_info ei;
		double weight = w;
		ei.weight = weight;
		gr.set_edge_info(e_new, ei);
		gr.set_edge_weight(e_new, weight);
		assert(e_new != NULL);
		assert(gr.ewrt.find(e_new) != gr.ewrt.end());
		assert(gr.einf.find(e_new) != gr.einf.end());
		assert(combinedRes.subpaths.size() >= 1);
		for(const path& p: combinedRes.subpaths)
		{
			assert(p.v.front() == source);
			assert(p.v.back() == target);
		}
	}
	assert(!gr.refine_splice_graph());

	return 0;
}

/* remove abutting edge, then call divide_conquer(source, target, res) again */
bool aster::divide_conquer_abutting(aster_index ai)
{	
	int s = ai.s();
	int t = ai.t();
	assert(s < gr.num_vertices() && t < gr.num_vertices() && s >= 0 && t >= 0);
	assert(s < t - 1);
	assert(gr.out_degree(s) >= 1 || gr.in_degree(t) >= 1);

	if(! gr.edge_exists(s,t)) return false;

	if(verbose >= 2) 
	{
		string msg = "aster processing subgraph, vertex [" + to_string(s) + ", " + to_string(t) + "]"; 
		// msg += " (topoIndex [" + to_string(source) + "," + to_string(target) + "]), ";
		msg += "removing the direct edge between them."; 
		cout << msg << endl;
	}

	// remove abutting edge & recusion
	edge_descriptor e = gr.edge(s, t).first;
	assert(e != null_edge);
	double w = gr.get_edge_weight(e);
	assert(gr.edge(e));
	aster_result res0 = edgeres.at(e);
	edgeres.erase(e);
	gr.remove_edge(e);
	assert(gr.check_path(s, t));
	divide_conquer(ai);

	// combine results
	PEB peb = gr.edge(s, t);
	assert(peb.second == true);
	edge_descriptor e1 = peb.first;
	double w1 = gr.get_edge_weight(e1);
	aster_result& res1 = edgeres.at(e1);
	aster_result rescomb;
	res_combine_parallel(res1, res0, rescomb);

	replace_aster_index_to_one_edge(ai, w + w1, rescomb);

	if(verbose >= 2) 
	{
		string msg = "aster processed subgraph, vertex [" + to_string(s) + ", " + to_string(t) + "]"; 
		// msg += " (topoIndex [" + to_string(source) + "," + to_string(target) + "]), ";
		msg += "with a direct edge between them."; 
		cout << msg << endl;
	}

	return true;
}	


bool aster::divide_conquer_cut_termini(aster_index ai)
{
	int s = ai.s();
	int t = ai.t();
	assert(s < gr.num_vertices() && t < gr.num_vertices() && s >= 0 && t >= 0);
	assert(s < t - 1);
	assert(gr.out_degree(s) >= 1 || gr.in_degree(t) >= 1);
	if(gr.edge_exists(s,t)) return false;
	if(gr.out_degree(s) == 1 || gr.in_degree(t) == 1) return false;
	assert(ai.size() >= 3);

	// get disjoint subgraphs' indices
	set<aster_index> subgraphIntervals;
	int subgraphNum = divide_conquer_cut_termini_find(ai, subgraphIntervals);
	if(subgraphNum <= 1 || subgraphIntervals.size() <= 1) return false;

	if(verbose >= 2) 
	{
		string msg = "aster processing subgraph, vertex [" + to_string(s) + ", " + to_string(t) + "] "; 
		msg += "splitting disjoint graphs at termini"; 
		msg += "\n\tsplitted subgraphIntervals:";
		cout << msg << endl;
		for(const auto& iv: subgraphIntervals)
		{
			string msg;
			msg += "\t";
			msg += "[" + to_string(iv.s()) + "-" + to_string(iv.t()) + "]: ";
			cout << msg;
			printv(iv.get_index());
			cout << endl;
		}
	}

	// D&C sub intervals INCLUDING s and t
	for(const aster_index& subItv: subgraphIntervals)
	{
		divide_conquer(subItv);
		// Now subItv is decomposed to one single edge
		// next iteration should remove the abutting edge before decomposing the new subItv
		for(int i : subItv.get_index()) assert(gr.degree(i) == 0 || i == s || i == t);
		assert(gr.edge_exists(s, t));
	}

	for(int i: ai.get_index()) assert(gr.degree(i) == 0 || i == s || i == t);
	assert(gr.edge_exists(s, t));
	// assert(gr.compute_num_paths(s, t, 2) == 1);
	assert(edgeres.find(gr.edge(s,t).first) != edgeres.end());

	if(verbose >= 2) 
	{
		string msg = "aster processed subgraph, vertex [" + to_string(s) + ", " + to_string(t) + "]"; 
		msg += "with " + to_string(subgraphIntervals.size()) +" disjoint graphs at termini"; 
		cout << msg << endl;
	}
	return true;
}

/* return size of interval pairs, INCLUDING original s and t
*  return -1 if not exists 
*  if source/target is cut vertex, 
* 		then all subgraph's sources must be source-vertex's out-edge targets
* 		then all subgraph's targets must be target-vertex's in-edge  sources
*/
int aster::divide_conquer_cut_termini_find(aster_index ai, set<aster_index>& aiSubIntervals)
{
	int s = ai.s();
	int t = ai.t();
	assert(s < gr.num_vertices() && t < gr.num_vertices() && s >= 0 && t >= 0);
	assert(s < t - 1);
	assert(gr.out_degree(s) >= 1 || gr.in_degree(t) >= 1);
	if(gr.edge_exists(s,t)) return -1;
	if(gr.out_degree(s) == 1 || gr.in_degree(t) == 1) return -1;
	if(ai.size() < 3) return -1;

	// build undirected graph w/o s and t & find connected components
	undirected_graph ug;
	ug.clear();
	map<int, int> ai2newi;
	for(int i = 1; i < ai.size() - 1; i++)
	{
		ug.add_vertex();
		ai2newi.insert({ai.at(i), i - 1});
	}

	for(int i = 1; i < ai.size() - 1; i++)
	{
		int k = ai.at(i);
		PEEI pei;
		pei = gr.in_edges(k);
		for(edge_iterator it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
		{
			edge_descriptor e = (*it1);
			int s = e->source();
			int t = e->target();
			if(! ai.index_found(s)) return -1;
			if(! ai.index_found(t)) return -1;
			if(s == ai.s() || t == ai.t()) continue;
			int news = ai2newi.at(s), newt = ai2newi.at(t);
			ug.add_edge(news, newt);
		}
		pei = gr.out_edges(k);
		for(edge_iterator it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
		{
			edge_descriptor e = (*it1);
			int s = e->source();
			int t = e->target();
			if(! ai.index_found(s)) return -1;
			if(! ai.index_found(t)) return -1;
			if(s == ai.s() || t == ai.t()) continue;
			int news = ai2newi.at(s), newt = ai2newi.at(t);
			ug.add_edge(news, newt);
		}
	}

	vector< set<int> > vv = ug.compute_connected_components();
	if(vv.size() <= 1) return -1;

	// insert cc to aiSubIntervals
	for(const auto& cc: vv)
	{
		vector<int> ccSort(cc.begin(), cc.end());
		for(int i = 0; i < ccSort.size(); i ++)
		{
			ccSort.at(i) = ai.at(ccSort.at(i) + 1);
		}
		sort(ccSort.begin(), ccSort.end());
		ccSort.insert(ccSort.begin(), ai.s());
		ccSort.push_back(ai.t());
		aiSubIntervals.insert(aster_index(ccSort));
	}

	// assertions and validations
	for(const auto& sub: aiSubIntervals)
	{
		assert(sub.size() >= 3);
		int subsource = sub.at(2);
		int subtarget = sub.at(sub.size() - 1 - 1);
		assert(subsource > s);
		assert(subtarget < t);
		for(int i: sub.get_index())
		{
			if(i == s || i == t) continue;
			PEEI peei;
			peei = gr.out_edges(i);
			for(edge_iterator it1 = peei.first, it2 = peei.second; it1 != it2; it1++)
			{
				assert(sub.index_found((*it1)->target()) || (*it1)->target() == t);
				assert(sub.index_found((*it1)->source()) || (*it1)->source() == s);
			}
			peei = gr.in_edges(i);
			for(edge_iterator it1 = peei.first, it2 = peei.second; it1 != it2; it1++)
			{
				assert(sub.index_found((*it1)->target()) || (*it1)->target() == t);
				assert(sub.index_found((*it1)->source()) || (*it1)->source() == s);
			}
		}
	}

	if(aiSubIntervals.size() == 0) return -1;
	return aiSubIntervals.size();
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

	if(verbose >= 2)
	{
		string msg = "aster processing subgraph, vertex [" + to_string(s) + ", " + to_string(t) + "]"; 
		msg += "with subgraphs at articulation point " + to_string(pivot);
		cout << msg << endl;
	}

	divide_conquer(aileft); 	
	PEB peb1 = gr.edge(s, pivot);
	edge_descriptor e1 = peb1.first;
	assert(peb1.second);

	divide_conquer(airight);
	PEB peb2 = gr.edge(pivot, t);
	edge_descriptor e2 = peb2.first;
	assert(peb2.second);

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
		msg += "with subgraphs at articulation point " + to_string(pivot);
		cout << msg << endl;
	}

	return true;
}

// find a pivot-TopoIndex s.t. removing this vertex will split grpah to two parts between [source, pivot] and [pivot, target]
// return -1 if cannot fine artivulation point
int aster::divide_conquer_articulation_find(aster_index ai, aster_index& left, aster_index& right) 
{
	int s = ai.s();
	int t = ai.t();
	if(ai.size() <= 2) return false;
	if(! gr.check_path(s, t)) return false;

	int pivot = -1;
	int index = -1;

	//TODO: optimize
	splice_graph grLocal(gr);
	for(int i = 0; i < grLocal.num_vertices(); i++)
	{
		if(ai.index_found(i)) continue;
		grLocal.clear_vertex(i);
	}

	for(int i = 1; i < ai.size() - 1; i++)
	{
		int v = ai.at(i);
		splice_graph gr2(grLocal);
		gr2.clear_vertex(v);
		if(gr2.check_path(s, t)) continue;
		pivot = v;
		index = i;
		break;
	}
	if (pivot >= t || pivot <= s || pivot <= 0) return -1;

	if(verbose >= 2) 
	{
		string msg = "\t found articulation point = " + to_string(pivot);
		cout << msg << endl;
	}

	// split ai
	ai.split(index, left, right);
	assert(left.t() == pivot);
	assert(right.s() == pivot);

	//assertion
	splice_graph gr2(gr);
	gr2.clear_vertex(pivot);
	for(int i = 0; i < gr2.num_vertices(); i++) if(!ai.index_found(i)) gr2.clear_vertex(i); 
	assert(! gr2.check_path(s, t)); 

	return pivot;
}

/*  combines subpaths of left and right sides of articulation point k
	res1, res2 could be empty
	If comb is non-empty, it is an addition
*/
bool aster::res_combine_consecutive(aster_result& res1,  aster_result& res2, aster_result& comb) const
{
	assert(res1.size() >= 1);
	assert(res2.size() >= 1);
	// if(res1.subpaths.size() == 0 && res2.subpaths.size() == 0) return false;
	// if(res1.subpaths.size() == 0) 
	// {
	// 	comb.subpaths.insert(comb.subpaths.end(), res2.subpaths.begin(), res2.subpaths.end());
	// 	comb.dist += res2.dist;
	// 	return true;
	// }
	// if(res2.subpaths.size() == 0) 
	// {
	// 	comb.dist += res1.dist;
	// 	comb.subpaths.insert(comb.subpaths.end(), res1.subpaths.begin(), res1.subpaths.end());
	// 	return true;
	// }

	assert(res1.subpaths.size() > 0);
	assert(res1.subpaths[0].v.size() > 0);
	assert(res2.subpaths.size() > 0);
	assert(res2.subpaths[0].v.size() > 0);
	
	int k = res1.subpaths[0].v.back();
	int combOriSize = comb.subpaths.size();

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

	int index1 = -1;
	int index2 = -1;
	if (strategy == aster_strategy::SHORT) 
	{
		index1 = find_shortest_path(res1);
		index2 = find_shortest_path(res2);
	} 
	else if (strategy == aster_strategy::LONG) 
	{
		index1 = find_longest_path(res1);
		index2 = find_longest_path(res2);
	} 
	else 
	{  // Default: strategy == aster_strategy::HEAVY)
		index1 = find_heaviest_path(res1);
		index2 = find_heaviest_path(res2);
	}
	
	assert(index1 >= 0);
	assert(index2 >= 0);
	assert(index1 < res1.subpaths.size());
	assert(index2 < res2.subpaths.size());
	
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

	assert(comb.subpaths.size() >= combOriSize + res1.subpaths.size() + res2.subpaths.size() - 1);

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

/*	combine res1 and res2, add them to comb.
	If comb is non-empty, it is an addition
*/ 
bool aster::res_combine_parallel(aster_result& res1,  aster_result& res2, aster_result& comb, bool frontSame, bool backSame) const
{	
	assert(res1.size() >= 1);
	assert(res2.size() >= 1);
	// if(res1.subpaths.size() == 0 && res2.subpaths.size() == 0) return false;
	// if(res1.subpaths.size() == 0) 
	// {
	// 	comb.subpaths.insert(comb.subpaths.end(), res2.subpaths.begin(), res2.subpaths.end());
	// 	comb.dist += res2.dist;
	// 	return true;
	// }
	// if(res2.subpaths.size() == 0) 
	// {
	// 	comb.subpaths.insert(comb.subpaths.end(), res1.subpaths.begin(), res1.subpaths.end());
	// 	comb.dist += res1.dist;
	// 	return true;
	// }
	assert(res1.subpaths.size() > 0);
	assert(res1.subpaths[0].v.size() > 0);
	assert(res2.subpaths.size() > 0);
	assert(res2.subpaths[0].v.size() > 0);

	int s = res1.subpaths[0].v.front();
	int t = res1.subpaths[0].v.back();

	for(const path& p: res1.subpaths)	
	{
		assert(p.v.size() >= 1); 
		if(frontSame) assert(p.v.front() == s);
		if(backSame)  assert(p.v.back() == t);
	}
	for(const path& p: res2.subpaths)	
	{
		assert(p.v.size() >= 1); 
		if(frontSame) assert(p.v.front() == s);
		if(backSame)  assert(p.v.back() == t);
	}
    
    // Collapse existing paths (if duplicated) by adding their counts
	map<vector<int>, double> unique_paths;
    for(const path& p : comb.subpaths) 
	{
        unique_paths[p.v] = p.abd;
    }
    for(const path& p : res1.subpaths) 
	{
        if(unique_paths.find(p.v) != unique_paths.end()) unique_paths[p.v] += p.abd;
        else unique_paths[p.v] = p.abd;
    }
    for(const path& p : res2.subpaths) 
	{
        if(unique_paths.find(p.v) != unique_paths.end()) unique_paths[p.v] += p.abd;  
        else unique_paths[p.v] = p.abd;
    }

    comb.subpaths.clear();
    for(const auto& pair : unique_paths) {
        comb.subpaths.push_back(path(pair.first, pair.second));
    }

	//TODO: calculate dist
    assert(comb.subpaths.size() >= 1);
    return true;
}

/* examine if dnc unitig between source to target; if true, populate res */
bool aster::divide_conquer_unitig(aster_index ai, splice_graph& localGr)
{
	int s = ai.s();
	int t = ai.t();
	assert(s < gr.num_vertices() && t < gr.num_vertices() && s >= 0 && t >= 0);
	assert(s < t);
	assert(gr.out_degree(s) >= 1 || gr.in_degree(t) >= 1);


	// local graph, num path >= 1
	assert(ai.size() >= 2);
	if(ai.size() ==  2) assert(gr.valid_path({s, t}) && localGr.valid_path({s, t}));
	int n = localGr.compute_num_paths(s, t, 2);
	assert(n >= 1);
	if(n > 1) return false;

	// find unitig in localGr
	int ss = s;
	vector<int> unitig;
	while(true)
	{
		unitig.push_back(ss);
		assert(ss <= t);
		if (ss == t) break;
		assert(gr.out_degree(ss) == 1 || ss == s || ss == t);
		assert(localGr.out_degree(ss) == 1);
		edge_descriptor e = (*(localGr.out_edges(ss).first));
		ss  = e->target();
	}
	assert(unitig.size() >= 2);
	assert(localGr.valid_path(unitig));

	// populate unitig
	aster_result unitigRes;
	bool   _avg_ = avgMode;       							// average if true, geom mean if false
	double w     = _avg_? 0: 1;
	double c = 0.0;
	{
		int ss = unitig.at(0);
		int tt = unitig.at(1);
		assert(ai.index_found(ss));
		assert(ai.index_found(tt));
		auto [e, edgeExists] = gr.edge(ss, tt);
		assert(edgeExists);
		double ew = gr.get_edge_weight(e);
		w = _avg_? (w + ew): (w * ew);
		c += 1.0;
		unitigRes = edgeres.at(e);
	}
	
	for(int i = 1; i < unitig.size() - 1; i++)
	{
		int ss = unitig.at(i);
		int tt = unitig.at(i + 1);
		assert(ai.index_found(ss));
		assert(ai.index_found(tt));
		auto [e, edgeExists] = gr.edge(ss, tt);
		assert(edgeExists);

		double ew = gr.get_edge_weight(e);
		w = _avg_? (w + ew): (w * ew);
		aster_result __res__;
		res_combine_consecutive(unitigRes, edgeres.at(e), __res__);
		unitigRes = __res__;
		c += 1.0;
	}
	w = _avg_? (w / c): pow(w, 1.0/ c);
	assert(valid_paths(unitigRes));

	if(verbose >= 2) 
	{
		string msg = "aster processed subgraph, vertex [" + to_string(s) + ", " + to_string(t) + "]"; 
		msg += "with a unitig path."; 
		cout << msg << endl;
		cout << "\t local path:";
		printv(unitig);
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

	// aster_result res;
	// double abd = gr.get_vertex_weight(s);
	throw runtime_error("single vertex should not be called");
	return false;
}

int aster::init_edgeres()
{
	// init edgeres for each edge
    edgeres.clear();

    PEEI pei = gr.edges();
    for(edge_iterator it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
    {
        edge_descriptor e = *it1;
        double w = gr.get_edge_weight(e);
        edgeres[e] = aster_result(vector<int>{e->source(), e->target()}, w);
    }
    
	// Then handle hyper-edges from hs.edges
    for(size_t i = 0; i < hs.edges.size(); i++) 
    {
        const vector<int>& hedge = hs.edges[i];
        if(hedge.size() < 2) continue;
        
        // Convert edge indices to node sequence
        vector<int> nodes;
        bool valid = true;
        
        // hyper_edge: Get first edge to get starting node
        edge_descriptor first_edge = i2e[hedge[0]];
        nodes.push_back(first_edge->source());
        nodes.push_back(first_edge->target());
        
        // hyper_edge:  Process remaining edges
        for(size_t j = 1; j < hedge.size(); j++) 
        {
            if(hedge[j] == -1) {
                valid = false;
                break;
            }
            
            edge_descriptor curr_edge = i2e[hedge[j]];
            
            // Verify connectivity
            if(curr_edge->source() != nodes.back()) {
                valid = false;
                break;
            }
            
            nodes.push_back(curr_edge->target());
        }
        if(!valid || nodes.size() < 2) continue;

		double w = hs.ecnts[i];
        
        // Add hyper-edge to graph: Check if edge already exists
		if(auto [e, exists] = gr.edge(nodes.front(), nodes.back()); exists)
		{
            // Edge exists - combine the aster_results
			assert (e != null_edge);
            aster_result new_res(nodes, w);
            aster_result combined_res;
            assert(edgeres.find(e) != edgeres.end());
            res_combine_parallel(edgeres[e], new_res, combined_res);
            edgeres[e] = combined_res;
            
            double old_w = gr.get_edge_weight(e);
            gr.set_edge_weight(e, old_w + w);  // Update edge weight
        } 
		else 
		{
            // Create new edge for this hyper-edge
            edge_descriptor new_edge = gr.add_edge(nodes.front(), nodes.back());
			edge_info ei;
			ei.weight = w;
			ei.strand = gr.strand;
			gr.set_edge_info(new_edge, ei);
            gr.set_edge_weight(new_edge, w);
            edgeres[new_edge] = aster_result(nodes, w);
        }
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
	// v2tp.clear();

	vector<bool> visited(gr.num_vertices(), false);
	for(int i = 0; i < gr.num_vertices(); i++)	
	{
		topological_sort_vertices_visit(i, visited, tp2v);
	}
	reverse(tp2v.begin(), tp2v.end());
	// v2tp.resize(tp2v.size());
	// for(int i = 0; i < tp2v.size(); i++) v2tp[tp2v[i]] = i;

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
/* 
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
*/

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

// clear empty vertex from graph
int aster::purge_empty_vertex()
{
	for(int i = 0; i < gr.num_vertices(); i++)
	{
		if(gr.get_vertex_info(i).type == EMPTY_VERTEX) gr.clear_vertex(i);
	}
	gr.refine_splice_graph();
	return 0;
}

int aster::remove_small_junctions()
{
	SE se;
	for(int i = 1; i < gr.num_vertices() - 1; i++)
	{
		if(gr.degree(i) <= 0) continue;

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

/* balance vertices, remove small junctions, refine graph, build index, build hs */
int aster::prepare_graph()
{
	for(int i = 1; i < gr.num_vertices() - 1; i++) balance_vertex(i);
	for(int i = 1; i < gr.num_vertices() - 1; i++) balance_vertex(i);
	for(int i = 1; i < gr.num_vertices() - 1; i++) balance_vertex(i);
	for(int i = 1; i < gr.num_vertices() - 1; i++) balance_vertex(i);
	// remove_small_junctions();
	purge_empty_vertex();
	gr.refine_splice_graph();
	gr.get_edge_indices(i2e, e2i);
	hs.build(gr, e2i);

	return 0;
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

// return index of heaviest path (highest abundance) in res.subpaths
int aster::find_heaviest_path(const aster_result& res) const
{
    int heaviestPathIndex = -1;
    double heaviestPathWeight = -1.0;
    if (res.subpaths.size() == 0) return heaviestPathIndex;

    heaviestPathIndex = 0;
    heaviestPathWeight = res.subpaths.front().abd;
    for(int i = 0; i < res.subpaths.size(); i++) 
    {
        const path& p = res.subpaths.at(i);
        if(p.abd <= heaviestPathWeight) continue;
        heaviestPathWeight = p.abd;
        heaviestPathIndex = i;
    }
    return heaviestPathIndex;
}

/* 
*  Assume internal nodes must be "closed" in [s, t]
*  only remove edges in [s,t] (toposorted). 
*  Maybe not necessarily removes all edges from/to vertices in [s,t] if it is not closed
*  i.e. for nodes reachable from s and reachable to t, they cannot reach any other nodes but themselves and s,t
*/
edge_descriptor aster::replace_aster_index_to_one_edge(aster_index ai, double w, aster_result res)
{
	assert_closed_vertex_interval(ai);
	int s = ai.s();
	int t = ai.t();
	assert(s < gr.num_vertices() && t < gr.num_vertices() && s >= 0 && t >= 0);
	assert(s < t);
	assert(gr.out_degree(s) >= 1 || gr.in_degree(t) >= 1);	
	assert(res.subpaths.size() >= 1);
	set<int> __set__;
	for(const path& p: res.subpaths)
	{
		assert(p.v.front() == s);
		assert(p.v.back() == t);
		for(int i: p.v) __set__.insert(i);
	}
	for(int i: __set__) assert(ai.index_found(i) || gr.degree(i) == 0);

	// under 'closed' assumption all edges in (s, t) open interval must have sources/targets in [s, t] closed interval
	for(int k: ai.get_index())
	{
		if(gr.degree(k) <= 0) continue;
		if(k == s || k == t) continue;
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
			assert(ai.index_found(e->source()));
			assert(ai.index_found(e->target()));
			edgeres.erase(e);
			gr.remove_edge(e);
		}
	}
	if(auto [e, b] = gr.edge(s,t); b) 
	{
		edgeres.erase(e);
		gr.remove_edge(e);
	}

	for(int i: ai.get_index()) assert(gr.degree(i) == 0 || i == ai.s() || i == ai.t());

	// put edge & res
	assert(!gr.edge_exists(s,t));
	edge_descriptor e_new = gr.add_edge(s, t);
	edge_info ei;
	ei.weight = w;
	gr.set_edge_info(e_new, ei);
	gr.set_edge_weight(e_new, w);
	assert(e_new != NULL);
	assert(gr.ewrt.find(e_new) != gr.ewrt.end());
	assert(gr.einf.find(e_new) != gr.einf.end());
	edgeres[e_new] = res;

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
	cout << "\t " << dnc_counter_resolve_trivial_node;
	cout << "\t " << dnc_counter_resolve_intersection_edge;
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
	// const vector<int>& v1 = p1.v; 
	// const vector<int>& v2 = p2.v;
	// assert(v1.size() > 0);
	// assert(v2.size() > 0);
	// if(v1.back() > v2.front()) return -1;
	// if(v2.back() > v1.front()) return -1;
	int edits;
	// edits = - basic_algo::ref_sw_query_nw(v1, v2, -1, -2, 0);
	return edits = -1;
}

// remove isolated vertices from aster_index, edit in-place
int aster::non_isolated_vertex_index(aster_index& ai) const
{
	while(true)
	{
		bool isErase = false;
		for(int i = 0; i < ai.size(); i++)
		{
			int v = ai.at(i);
			if(gr.degree(v) >= 1) continue;
			if(i == 0 || i == ai.size() - 1) 
			{
				throw runtime_error("subgraph source/sink should not be isolated");
			}
			ai.erase_index(i);
			isErase = true;
			break;
		}
		if(isErase) continue;
		else break;
	}
	return 0;
}

// all internal nodes (excl. s and t) must have internal edges whose target/source can be found in aster_index
int aster::assert_closed_vertex_interval(const aster_index &ai)
{
	if(ai.size() <= 2) return 0;
	int s = ai.s();
	int t = ai.t();
	for(int i = 1; i < ai.size() - 1; i++)
	{
		if (gr.degree(i) <= 0) continue;
		int k = ai.at(i);

		PEEI pi = gr.in_edges(k);
		for(edge_iterator it = pi.first; it != pi.second; it++)
		{
			edge_descriptor e = *it;
			assert(gr.edge(e));
			assert(ai.index_found(e->source()) || e->target() == t);
		}
		PEEI po = gr.out_edges(k);
		for(edge_iterator it = po.first; it != po.second; it++)
		{
			edge_descriptor e = *it;
			assert(gr.edge(e));
			assert(ai.index_found(e->target()) || e->source() == s);
		}
	}
	return 0;
}
