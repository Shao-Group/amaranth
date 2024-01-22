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

	divide_conquer();
	return 0;
}

int aster::divide_conquer()
{
	assert(gr.num_vertices() > 2);
	int s = 0;
	int t = gr.num_vertices() - 1;
	aster_result res;
	divide_conquer(s, t, res);
	paths = res.subpaths;
	assert(paths.size() > 0);
	for(const path & p : paths)	assert(gr.valid_path(p.v));
	return 0;
}

// divide_conquer(i ,j) solves a subproblem between 
int aster::divide_conquer(int source, int target, aster_result& res)
{
	assert(source <= target);
	int s = source;
	int t = target;
	assert(res.subpaths.size() == 0);
	assert(res.dist == -1);

	if (divide_conquer_single_vertex(s, t, res))      return 0;
	if (divide_conquer_unitig(s, t, res))             return 0;
	if (divide_conquer_disjoint(s, t, res))     	  return 0;
	if (divide_conquer_abutting(s, t, res))			  return 0;
	
	assert(0);
	return -1;
}

/* examine if dnc unitig between s to t; if true, populate res */
bool aster::divide_conquer_unitig(int s, int t, aster_result& res)
{
	assert(s < t);
	assert(gr.out_degree(s) >= 1 || gr.in_degree(t) >= 1);

	if(gr.out_degree(s) > 1 || gr.in_degree(t) > 1) return false;

	int n = gr.compute_num_paths(s, t, 2);
	assert(n >= 1);
	if(n > 1) return false;

	vector<int> vertexPath;
	gr.compute_shortest_path(s, t, vertexPath);
	assert(vertexPath.size() > 0);
	double abd = 0;															//TODO
	res.subpaths.push_back(path(vertexPath, abd));
	res.dist = 0;
	
	return true;
}

/* examine if dnc single vertex; if true, populate res */
bool aster::divide_conquer_single_vertex(int s, int t, aster_result& res)
{
	assert(s <= t);
	if(s != t) return false;
	res.subpaths.push_back({}); 
	res.dist = 0;
	return true;
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
		if (i < gr.num_vertices() - 1) 		// proceeding out edges
		{
			PEEI ei = gr.out_edges(i);
			edge_iterator it1 = ei.first; 
			edge_iterator it2 = ei.second;
			for(; it1 != it2; it1++) assert(gr.get_vertex_info(i).rpos <= gr.get_vertex_info((*it1)->target()).lpos);
		}
		if (i >= 1) 						// trailing in edges
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
	gr.refine_splice_graph();
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


int aster::get_transcripts()
{
	//TODO:
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


// exponential penalty guarantees to violate triangle inequality
int aster::event_size_penalty(int eventSize)
{
	assert(eventSize >= 0);
	return pow(2, eventSize) - 1;
}

/*
* return: 
*		positive int: number of edits
*		assertion error: size not positive
*/
int aster::path_distance(const path& p1, const path& p2)
{
	const vector<int>& v1 = p1.v; 
	const vector<int>& v2 = p2.v;
	assert(v1.size() > 0);
	assert(v2.size() > 0);
	if(v1.back() > v2.front()) return -1;
	if(v2.back() > v1.front()) return -1;

	int edits = - basic_algo::ref_sw_query_nw(v1, v2, -1, -2, 0);

	return edits; //TODO: not finished
}

