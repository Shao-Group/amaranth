/*
Part of Aster Assembler
(c) 2024 by Xiaofei Carl Zang, Mingfu Shao, and The Pennsylvania State University.
See LICENSE for licensing.
*/

#ifndef __ASTER_H__
#define __ASTER_H__

#include "splice_graph.h"
#include "hyper_set.h"
#include "equation.h"
#include "router.h"
#include "path.h"
#include "config.h"

typedef map< edge_descriptor, vector<int> > MEV;
typedef pair< edge_descriptor, vector<int> > PEV;
typedef pair< vector<int>, vector<int> > PVV;
typedef pair<PEE, int> PPEEI;
typedef map<PEE, int> MPEEI;
typedef pair<int, int> PI;
typedef map<int, int> MI;
// typedef vector<path> aster_dp_dot;
// typedef vector<aster_dp_dot> aster_dp_row;
// typedef vector<aster_dp_row> aster_dp_table; 

struct aster_result
{
	vector<path> subpaths;
	int dist;
};

/* aster class performs the preparation and revise work */
class aster
{
public:
	aster(const splice_graph &gr, const hyper_set &hs);
	int assemble();
	static int print_stats();

// static stats
private:
	inline static int num_graph = 0;
	inline static int num_intersecting_graph = 0;
	inline static int num_intersecting_graph2 = 0;
	inline static int num_intron = 0;
	inline static int num_exon = 0;
	inline static int num_intersecting_intron_count = 0;
	inline static int num_intersecting_intron_pair = 0;

public:
	splice_graph gr;					// splice graph
	hyper_set hs;						// hyper edges
	MEI e2i;							// edge map, from edge to index, sorted
	VE i2e;								// edge map, from index to edge, sorted
	vector<path> paths;					// predicted paths
    vector<transcript> trsts;			// predicted transcripts
	vector<transcript> non_full_trsts;		// predicted non full length transcripts

private: 
	int topological_sort_vertices();
	int topological_sort_index_edges();
	int aggressive_purge_intersecting_edges();
	int balance_vertex(int);

	int divide_conquer();
	int divide_conquer(int source, int target, aster_result& table);
	
	int get_transcripts();
	int make_stats();
};

/* astron class executes algorithms under clear assumptions */
class astron
{
public:
	astron(aster*, const vector<int>& _canons, const vector<int>& _illegal, const vector<int>& _alts = {}, string algo = "dp");

private:
	aster* as;
    vector<int> canons;			// canonical events
	vector<int> illegals;		// illegal events
    vector<int> alternatives;	// alternative events
	string aster_algo;

public:
	int dist;
	vector<path> paths;

private:
	int classify();
	int dynamic_programming();												// algo dp	//TODO:
	int divide_and_conquer();												// algo dnc //TODO:
	int dnc_combine(const vector<path> subpaths, int eventOfConcern);		// algo dnc //TODO:
	int heuristic();														// algo heuristic	//TODO:
	bool heuristic_search(vector<vector<int>>& ppNodes, int maxDist);
	bool heuristic_search(vector<edge_descriptor>& edges, int maxDist);
	int greedy();
	bool greedy_longest_path();
	bool greedy_edit_path(int maxDist);

	int collect_trivial_path();
	int event_size_penalty(int eventSize);
	int closest_path(vector<int> nodes, int maxDist);						// find closest path not exceeding max penalty
	int path_distance(const vector<int>& v1, const vector<int>& v2);
};

#endif