/*
Part of Aster Transcript Assembler
(c) 2024 by Xiaofei Carl Zang, Mingfu Shao, and The Pennsylvania State University.
See LICENSE for licensing.
*/

#include <algorithm>
#include <cmath>
#include "aster.h"
#include "basic_algo.h"

aster::aster(const splice_graph &g, const hyper_set &h, bool r)
	: gr(g), hs(h), random_ordering(r)
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
}

int aster::assemble()
{
	vector<int> canonical;
	vector<int> illegal;
	canonical.clear();
	illegal.clear();
	astron asterPetal(this, canonical, illegal);
	return 0;
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

int aster::make_stats()
{
	// num graph, exon, intron
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
	
	// num intersecting introns, intron_pairs, graphs
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
astron::astron(const aster* _as, const vector<int>& _canons,  const vector<int>& _illegal, const vector<int>& _alts, string _algo)
	: as(_as), canons(_canons), illegals(_illegal), alternatives(_alts), aster_algo(_algo), dist(-1)
{
	classify();

	aster_algo = "heuristic";
	if 		(aster_algo == "dp") 			dist = dynamic_programming();
	else if (aster_algo == "dnc")  			dist = divide_and_conquer();
	else if (aster_algo == "heuristic") 	dist = heuristic();
	else assert(0);


	// if(alternatives.size() == 0) collect_trivial_path();
	// assert(dist >= 0);
	// assert(paths.size() >= 1);
}

int astron::classify()
{	
	for (int i: canons)    assert(find(illegals.begin(),  illegals.end(), i)  == illegals.end());
	for (int i: illegals)  assert(find(canons.begin(),    canons.end(), i)    == canons.end());
	if(alternatives.size() != 0) return 0;
	
	int altSize = as->gr.num_vertices() - 2 - canons.size() - illegals.size();
	assert(altSize >= 0);
	alternatives.resize(altSize);
	for (int i = 1; i < as->gr.num_vertices() - 1; i++)
	{
		if(find(illegals.begin(),  illegals.end(), i)  != illegals.end()) continue;
		if(find(canons.begin(),    canons.end(), i)    != canons.end())  continue;
		alternatives.push_back(i);
	}	
	return 0;
}

int astron::divide_and_conquer()
{
	int minDist = event_size_penalty(alternatives.size());
	
	for(int eventOfConcern: alternatives)
	{	
		int exclusiveEvents = 1;
		int eventPenalty = event_size_penalty(exclusiveEvents);
		unique_ptr<astron> canonChild; // TODO: one child is enough. Alt child distance can be calculated from canonChild

		assert(canonChild->dist >= 1);
		// dnc_combine(); // TODO:
		int newDist = canonChild->dist + eventPenalty;
		if (minDist > newDist) minDist = newDist;
	}

	return minDist;
}

int astron::dnc_combine(const vector<path> subpaths, int eventOfConcern)
{

}


/*
* Search for a new path using a new phasing path
* s.t. 1. the phasing path is completely contained
*	   2. the new path has edit distance < `x` from at least 1 of acceptable path
*	   	2.1 for (x = 1 to 5) if new search found, reset x, otherwise increase x
*/
int astron::heuristic()
{
	// sort pashing paths by counts
	vector<vector<int>> ppNodes;
	vector<int> 		ppCounts;
	as->hs.sort_nodes(ppNodes, ppCounts);
	assert(ppNodes.size()  == as->hs.nodes.size());
	assert(ppCounts.size() == as->hs.nodes.size());
	
	vector<edge_descriptor> edges;
	for(const auto& edgeIndex: as->e2i) edges.push_back(edgeIndex.first);

	// heuristics to cover a pp or edge through
	while (ppNodes.size() != 0 && edges.size() != 0)  
	{
		if (heuristic_search(ppNodes)) continue;
		else heuristic_search(edges);
	}
	

	return dist;
}

/* try to cover a pp through heuristics
*/
bool astron::heuristic_search(vector<vector<int>>& ppNodes)
{

}

/* try to cover an edge through heuristics
*/
bool astron::heuristic_search(vector<edge_descriptor>& edges)
{

}

int astron::dynamic_programming()
{

}

/*
* collect paths based on canon events
*/
int astron::collect_trivial_path()
{
	const splice_graph& gr = as->gr;

	vector<int> v = canons;
	sort(v.begin(), v.end());

	// v does not contain source & sink
	assert(v[0] != 0); 
	assert(v[v.size() - 1] != gr.num_vertices()); 
	assert(gr.valid_path(v)); //FIXME: not implemented
	

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

int astron::closest_path(vector<int> nodes)
{
	int x = 1; //TODO:
	for(int i = 0; i < paths.size(); i++)
	{
		int d = path_distance(paths[i].v, nodes, x);
		assert(d >= 0);
	}
}

/*
* return: 
*		positive int: number of edits
*		assertion error: size not positive
*/
int astron::path_distance(const vector<int>& v1, const vector<int>& v2, int maxAllow)
{
	assert(v1.size() > 0);
	assert(v2.size() > 0);
	if(v1.back() > v2.front()) return -1;
	if(v2.back() > v1.front()) return -1;

	int edits = 0;
	int k1 = v1.size() - 1;
	int k2 = v2.size() - 1;
	while(k1 >= 0 && k2 >= 0 && k1 <= v1.size() - 1 && k2 <= v2.size() - 1)
	{
		basic_algo::smith_waterman //FIXME:
	}


	return -1;
}



