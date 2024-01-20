/*
Part of Aster Transcript Assembler
(c) 2024 by Xiaofei Carl Zang, Mingfu Shao, and The Pennsylvania State University.
See LICENSE for licensing.
*/

#include <algorithm>
#include <cmath>
#include "aster.h"
#include "basic_algo.h"

aster::aster(const splice_graph &g, const hyper_set &h)
	: gr(g), hs(h)
{
	// prepare
    topological_sort_vertices();
	topological_sort_index_edges();
	make_stats();
	
	// revise
	aggressive_purge_intersecting_edges();
	topological_sort_vertices();
	topological_sort_index_edges();

	assemble();
	get_transcripts();
}

int aster::assemble()
{	

	if(gr.num_edges() == 0) return 0;
	if(gr.num_edges() == 2) return 0;
	assert(gr.num_vertices() > 2);

	//CLEAN: balance?
	if (true)
	{
		for(int i = 1; i < gr.num_vertices() - 1; i++) balance_vertex(i);
		for(int i = 1; i < gr.num_vertices() - 1; i++) balance_vertex(i);
	}

	dynamic_programming();
	return 0;
}

{
	int m = gr.num_vertices();
	aster_dp_table optPaths(m, aster_dp_row(m)); 	// opt[i][j] means optimal min evo dist aster_dp_dot from vertex i to j
	
	{
			else	    dynamic_programming(j, i, optPaths);
		}
	}
	return 0;
}

int aster::dynamic_programming(int source, int target, aster_dp_table& optPaths)
{
	assert(source <= target);
	int s = source;
	int t = target;

	assert(optPaths[s][t].size() == 0);
	if(s == t) optPaths[s][t].push_back({});

	if(int n = gr.compute_num_paths(s, t, 2); n <= 1) 
	{
		
		vector<int> vertexPath;
		gr.compute_shortest_path(s, t, vertexPath);
		assert(vertexPath.size() > 0);
		double abd = 0;
		paths.push_back(path(vertexPath, abd));
		return 0;
	}

}

/*
* sort vertices according to pair<lpos, rpos>
* assuming no two vertices have the same <lpos, rpos>
*/
int aster::topological_sort_vertices()
{
	assert(gr.num_vertices() >= 2);

	for(int i = 0; i < gr.num_vertices(); i++)	
	{
		if(i >= 1) assert(gr.get_vertex_info(i - 1).rpos <= gr.get_vertex_info(i).lpos);
		
		// proceeding out edges
		if (i < gr.num_vertices() - 1) 
		{
			PEEI ei = gr.out_edges(i);
			edge_iterator it1 = ei.first; 
			edge_iterator it2 = ei.second;
			for(; it1 != it2; it1++) assert(gr.get_vertex_info(i).rpos <= gr.get_vertex_info((*it1)->target()).lpos);
		}

		// trailing in edges
		if (i >= 1) 
		{
			PEEI ei = gr.in_edges(i);
			edge_iterator it1 = ei.first; 
			edge_iterator it2 = ei.second;
			for(; it1 != it2; it1++) assert(gr.get_vertex_info(i).lpos >= gr.get_vertex_info((*it1)->source()).rpos);
		}
	}

	return 0;
}

/*
* For future compatibility, edges should be sorted & fetched independent of splice_graph implementation
*/ 
int aster::topological_sort_index_edges()
{
	map<pair<int, int>, edge_descriptor> sortedEdges;
	PEEI pei = gr.edges(); 
	for(edge_iterator it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
	{
		edge_descriptor e = *it1;
		sortedEdges.insert({{e->source(), e->target()}, e});	
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
int aster::aggressive_purge_intersecting_edges()
{
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
			break;
		}
	}
	return 0;
}

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

int aster::make_stats()
{
	num_graph ++;
	for(int i = 1; i < gr.num_vertices() - 1; i++)
	{
		if(gr.degree(i) == 0) continue;
		num_exon ++;
	}
	for (int i = 0; i < i2e.size(); i ++) 
	{
		if(i2e[i]->source() == 0) continue; 
		if(i2e[i]->target() == gr.num_vertices() - 1) continue;
		num_intron++;
	}
	
	// num intersecting: introns, intron_pairs, graphs
	bool intersecting = false;
	for (int i = 0; i < i2e.size() - 1; i ++)
	{
		bool intersecting_edge = false;
		if(i2e[i] == null_edge) continue;
		for (int j = i + 1; j < i2e.size(); j ++)
		{
			if (i2e[j] == null_edge) continue;
			if (!gr.intersect(i2e[i], i2e[j])) continue;
			num_intersecting_intron_pair ++;
			intersecting_edge = true;
			intersecting = true;
		}
		
		if (!intersecting_edge) continue;
		num_intersecting_intron_count++;
		if(i == i2e.size() - 2) num_intersecting_intron_count ++;
	}

	if(intersecting) num_intersecting_graph ++;
	if(!gr.check_nested()) num_intersecting_graph2 ++;

	if(num_graph % 100 == 0)	print_stats();
	return 0;
}

int aster::print_stats()
{
	cout << "aster print stats" << endl;
	cout << "\t num graph " << num_graph  << endl; 
	cout << "\t num_intersecting_graph " << num_intersecting_graph  << endl; 
	cout << "\t num_intersecting_graph2 " << num_intersecting_graph2  << endl; 
	cout << "\t num intron " << num_intron << endl;
	cout << "\t num_exon " << num_exon  << endl; 
	cout << "\t num_intersecting_intron_count " << num_intersecting_intron_count << endl;
	cout << "\t num_intersecting_intron_pair " << num_intersecting_intron_pair << endl;
	cout << "aster printed stats" << endl;
	return 0;
}


/*
* astron is the divide-and-conquer base of aster
* aston runs an iterative choice of event of concern and remove it from the graph
*/
astron::astron(aster* _as, const vector<int>& _canons,  const vector<int>& _illegal, const vector<int>& _alts, string _algo)
	: as(_as), canons(_canons), illegals(_illegal), alternatives(_alts), aster_algo(_algo), dist(-1)
{
	classify();
	if 		(aster_algo == "dp") 			dist = dynamic_programming();
	// else if (aster_algo == "dnc")  			dist = divide_and_conquer();
	// else if (aster_algo == "heuristic") 	dist = heuristic();
	// else if (aster_algo == "greedy")		dist = greedy();
	else assert(0);

	// if(alternatives.size() == 0) collect_trivial_path();
	// assert(dist >= 0);
	// assert(paths.size() >= 1);
}

// int astron::classify()
// {	
// 	for (int i: canons)    assert(find(illegals.begin(),  illegals.end(), i)  == illegals.end());
// 	for (int i: illegals)  assert(find(canons.begin(),    canons.end(), i)    == canons.end());
// 	if(alternatives.size() != 0) return 0;
	
// 	int altSize = as->gr.num_vertices() - 2 - canons.size() - illegals.size();
// 	assert(altSize >= 0);
// 	alternatives.resize(altSize);
// 	for (int i = 1; i < as->gr.num_vertices() - 1; i++)
// 	{
// 		if(find(illegals.begin(),  illegals.end(), i)  != illegals.end()) continue;
// 		if(find(canons.begin(),    canons.end(), i)    != canons.end())  continue;
// 		alternatives.push_back(i);
// 	}	
// 	return 0;
// }

// int astron::divide_and_conquer()
// {
// 	int minDist = event_size_penalty(alternatives.size());
	
// 	for(int eventOfConcern: alternatives)
// 	{	
// 		int exclusiveEvents = 1;
// 		int eventPenalty = event_size_penalty(exclusiveEvents);
// 		unique_ptr<astron> canonChild; // TODO: one child is enough. Alt child distance can be calculated from canonChild

// 		assert(canonChild->dist >= 1);
// 		// dnc_combine(); // TODO:
// 		int newDist = canonChild->dist + eventPenalty;
// 		if (minDist > newDist) minDist = newDist;
// 	}

// 	return minDist;
// }

// int astron::dnc_combine(const vector<path> subpaths, int eventOfConcern)
// {

// }

// int astron::greedy()
// {	
// 	if(as->gr.num_edges() == 0) return 0;
// 	for(int i = 1; i < as->gr.num_vertices() - 1; i++) as->balance_vertex(i);
// 	for(int i = 1; i < as->gr.num_vertices() - 1; i++) as->balance_vertex(i);

// 	// int cnt = 0;
// 	// int n1 = paths.size();
// 	// while(true)
// 	// {
// 	// 	VE v;
// 	// 	double w = gr.compute_maximum_path_w(v);
// 	// 	if(w <= min_transcript_coverage) break;

// 	// 	int e = split_merge_path(v, w);
// 	// 	collect_path(e);
// 	// 	cnt++;
// 	// }
// 	// int n2 = paths.size();
// 	// if(verbose >= 2) printf("greedy decomposing produces %d / %d paths\n", n2 - n1, n2);
	
// 	greedy_longest_path();

// 	constexpr int maxDistAllow = 5;
// 	for(int i = 0; i < 5; i++)
// 	{
// 		if (greedy_edit_path(i)) continue;
// 		else break;
// 	}

// 	return 0;
// }

// bool astron::greedy_longest_path()
// {
// 	const splice_graph& gr = as->gr;

// 	for(int i = 0; i < gr.num_vertices(); i++) {} 


// }


// /*
// * Search for a new path using a new phasing path
// * s.t. 1. the phasing path is completely contained
// *	   2. the new path has edit distance < `x` from at least 1 of acceptable path
// *	   	2.1 for (x = 1 to 5) if new search found, reset x, otherwise increase x
// */
// int astron::heuristic()
// {
// 	// sort pashing paths by counts
// 	vector<vector<int>> ppNodes;
// 	vector<int> 		ppCounts;
// 	as->hs.sort_nodes(ppNodes, ppCounts);
// 	assert(ppNodes.size()  == as->hs.nodes.size());
// 	assert(ppCounts.size() == as->hs.nodes.size());
	
// 	vector<edge_descriptor> edges;
// 	for(const auto& edgeIndex: as->e2i) edges.push_back(edgeIndex.first);

// 	// heuristics to cover a pp or edge through
// 	while (ppNodes.size() != 0 && edges.size() != 0)  
// 	{
// 		if (heuristic_search(ppNodes)) continue;
// 		else if (heuristic_search(edges)) continue;
// 		break;
// 	}
	

// 	return dist;
// }

// /* try to cover a pp through heuristics
// */
// bool astron::heuristic_search(vector<vector<int>>& ppNodes, int maxDist)
// {
// 	for(auto it = ppNodes.begin(); it != ppNodes.end(); it ++)
// 	{
// 		const vector<int>& pp = *it;
// 		int pathIndex = closest_path(pp, 1);
// 		if (pathIndex < 0) continue;
// 		assert(pathIndex >= 0 && pathIndex < paths.size());
// 	}
// }

// /* try to cover an edge through heuristics
// */
// bool astron::heuristic_search(vector<edge_descriptor>& edges, int maxDist)
// {

// }

int astron::dynamic_programming()
{

}

/* collect paths based on canon events */
int astron::collect_trivial_path()
{
	const splice_graph& gr = as->gr;

	vector<int> v = canons;
	sort(v.begin(), v.end());

	// v does not contain source & sink
	assert(v[0] != 0); 
	assert(v[v.size() - 1] != gr.num_vertices()); 
	assert(gr.valid_path(v));

	// filter empty-vertex
	bool empty = false;
	for(int i = 0; i < v.size(); i++)
	{
		if(as->gr.get_vertex_info(v[i]).type != EMPTY_VERTEX) continue;
		empty = true;
		break;
	}

	path p;
	p.v = v;
	p.nf = empty? 1: 0;
	// p.abd = gr.get_edge_weight(i2e[e]);//FIXME:
	paths.push_back(p);

	return 0;
}

int astron::event_size_penalty(int eventSize)
{
	if (eventSize == 0) return 0;

	assert(eventSize >= 1);
	return pow(2, eventSize - 1);
}

// find closest path not exceeding maxDist
int astron::closest_path(vector<int> nodes, int maxDist) // FIXME:
{
	for(int i = 0; i < paths.size(); i++)
	{
		int d = path_distance(paths[i].v, nodes);
		assert(d >= 0);
	}
}

/*
* return: 
*		positive int: number of edits
*		assertion error: size not positive
*/
int astron::path_distance(const vector<int>& v1, const vector<int>& v2)
{
	assert(v1.size() > 0);
	assert(v2.size() > 0);
	if(v1.back() > v2.front()) return -1;
	if(v2.back() > v1.front()) return -1;

	int edits = - basic_algo::ref_sw_query_nw(v1, v2, -1, -2, 0);

	return edits; //TODO: not finished
}

