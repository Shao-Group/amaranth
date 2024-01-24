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
typedef vector<int> VI;
enum class comb_strat {GREEDY_MIN, GREEDY_MAX};

struct aster_result
{
	vector<path> subpaths;	// predicted paths, original v index, inclusive
	int dist = -1;
};

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
	inline static int dnc_counter_single = 0;
	inline static int dnc_counter_unitig = 0;
	inline static int dnc_counter_abutting = 0;
	inline static int dnc_counter_nested = 0;
	inline static int dnc_counter_disjoint = 0;

	const splice_graph& origr;			// original splice graph
	splice_graph gr;					// splice graph with modification
	hyper_set hs;						// hyper edges
	VI tp2v;							// DFS-based topologically sorted index to vertex index. This guarantees all disjoint subgraphs are gathered together
	VI v2tp;							// vertex index to DFS-based topologically sorted index. This guarantees all disjoint subgraphs are gathered together
	aster_mode mode;


public:
	vector<path> paths;					// predicted paths, original v index, inclusive
    vector<transcript> trsts;			// predicted transcripts, original v index
	vector<transcript> non_full_trsts;		// predicted non full length transcripts

private: 
	int topological_sort_vertices();
	int topological_sort_vertices_visit(int i, vector<bool>& visited);
	int topological_sort_index_edges();
	// bool aggressive_purge_intersecting_edges();
	int balance_vertex(int);

	int divide_conquer();
	int divide_conquer(int source, int target, aster_result& res);
	bool divide_conquer_single_vertex(int source, int target, aster_result& res);
	bool divide_conquer_unitig(int source, int target, aster_result& res);
	bool divide_conquer_abutting(int source, int target, aster_result& res);
	bool divide_conquer_disjoint_at_termini(int source, int target, aster_result& res);
	bool divide_conquer_disjoint_at_pivot(int source, int target, aster_result& res, comb_strat st);
	int  divide_conquer_find_pivot(int source, int target);
	int  divide_conquer_combine(aster_result& r1,  aster_result& r2, int pivot, aster_result& comb, comb_strat st);
	bool resolve_trivial_intersection(int source, int target, aster_result& res);

	int event_size_penalty(int eventSize);
	int path_distance(const path& p1, const path& p2);
	int edge_path_to_vertex_path(const VE& edgePath, VI& vertexPath);
	int find_longest_path(const aster_result& res);
	int find_shortest_path(const aster_result& res);
	int get_transcripts();
	
	string tp2v_to_string();
	int make_stats();
};

#endif