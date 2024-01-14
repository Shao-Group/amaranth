/*
Part of Aster Assembler
(c) 2024 by Xiaofei Carl Zang, Mingfu Shao, and The Pennsylvania State University.
See LICENSE for licensing.
*/

#include <algorithm>
#include <cmath>
#include "aster.h"

aster::aster(const splice_graph &g, const hyper_set &h, bool r)
	: gr(g), hs(h), random_ordering(r)
{
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
** sort vertices according to pair<lpos, rpos>
** assuming no two vertices have the same <lpos, rpos>
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
			PEEI ei = gr.in_edges(i);
			edge_iterator it1 = ei.first; 
			edge_iterator it2 = ei.second;
			for(; it1 != it2; it1++) assert(gr.get_vertex_info(i).rpos <= gr.get_vertex_info((*it1)->target()).lpos);
		}

		// trailing in edges
		if (i >= 1) 
		{
			PEEI ei = gr.out_edges(i);
			edge_iterator it1 = ei.first; 
			edge_iterator it2 = ei.second;
			for(; it1 != it2; it1++) assert(gr.get_vertex_info(i).lpos >= gr.get_vertex_info((*it1)->target()).rpos);
		}
	}

	return 0;
}

/*
** For future compatibility, edges should be sorted & fetched independent of splice_graph implementation
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
	assert(e2i.size() == index + 1 && i2e.size() == index + 1);

	return 0;
}


/*
** astron is the divide-and-conquer base of aster
** aston runs an iterative choice of event of concern and remove it from the graph
*/
astron::astron(const aster* _as, const vector<int>& _canonial,  const vector<int>& _illegal)
	: as(_as), canonial(_canonial), illegal(_illegal), dist(-1)
{
	classify();
	dist = divide_and_conquer();
	assert(dist >= 0);
}

int astron::classify()
{	
	for (int i: canonial) assert(find(illegal.begin(),  illegal.end(), i)  == illegal.end());
	for (int i: illegal)  assert(find(canonial.begin(), canonial.end(), i) == canonial.end());
	for (int i = 0; i< as->gr.num_vertices(); i++)
	{

	}	
}

int astron::divide_and_conquer()
{
	int minDist = event_size_penalty(alternative.size());

	for(int eventOfConcern: alternative)
	{	
		int exclusiveEvents = 1;
		int eventPenalty = event_size_penalty(exclusiveEvents);
		unique_ptr<astron> canonChild; // TODO: one child is enough. This child do work by removing all events of concern
		unique_ptr<astron> illegalChild;
		assert(canonChild->dist >= 1);
		assert(illegalChild->dist >= 1);
		// dnc_combine(); // TODO:
		int newDist = canonChild->dist + illegalChild->dist + eventPenalty;
		if (minDist > newDist) minDist = newDist;
	}

	return minDist;
}

int astron::dnc_combine(const vector<path> subpaths, int eventOfConcern)
{

}

int astron::event_size_penalty(int eventSize)
{
	assert(eventSize >= 1);
	return pow(2, eventSize - 1);
}

