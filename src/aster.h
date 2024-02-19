/*
Part of Aster Assembler
(c) 2024 by Xiaofei Carl Zang, Mingfu Shao, and The Pennsylvania State University.
See LICENSE for licensing.
*/

#ifndef __ASTER_H__
#define __ASTER_H__

#include <stdexcept>
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

struct aster_result
{
	vector<path> subpaths;	// predicted paths, original v index, inclusive
	int dist = -1;
	inline void clear() {subpaths.clear(); dist = -1;}
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
	inline static int num_intron = 0;
	inline static int num_exon = 0;
	inline static int num_intersecting_intron_count = 0;
	inline static int num_intersecting_intron_pair = 0;
	inline static int dnc_counter_single = 0;
	inline static int dnc_counter_unitig = 0;
	inline static int dnc_counter_abutting = 0;
	inline static int dnc_counter_cut_vertex = 0;
	inline static int dnc_counter_articulation_point_disjoint = 0;
	inline static int dnc_counter_resolve_trivial_intersection = 0;
	inline static int dnc_counter_resolve_trivial_paths = 0;
	inline static int dnc_counter_resolve_intersection_edge = 0;

	

	const splice_graph& origr;			// original splice graph
	splice_graph gr;					// splice graph with modification
	hyper_set hs;						// hyper edges
	VI tp2v;							// DFS-based topologically sorted index to vertex index. This guarantees all disjoint subgraphs are gathered together
	VI v2tp;							// vertex index to DFS-based topologically sorted index. This guarantees all disjoint subgraphs are gathered together
	aster_mode mode;

public:
	vector<path> paths;							// predicted paths, original v index, inclusive
    vector<transcript> trsts;					// predicted transcripts, original v index
	vector<transcript> non_full_trsts;			// predicted non full length transcripts
	bool successStatus = false;					// whether splice graph is decomposed successfully

private: 
	int topological_sort_vertices();
	int topological_sort_vertices_visit(int i, vector<bool>& visited, vector<int>& sorted);
	int topological_sort_index_edges();
//	bool aggressive_purge_intersecting_edges();
	int balance_vertex(int);
	int remove_small_junctions();

	int divide_conquer();
	int divide_conquer(int source, int target, aster_result& res);
	bool divide_conquer_single_vertex(int source, int target, aster_result& res);
	bool divide_conquer_unitig(int source, int target, aster_result& res);
	bool divide_conquer_abutting(int source, int target, aster_result& res);
	bool divide_conquer_cut_termini(int source, int target, aster_result& res);
	int  divide_conquer_cut_termini_find(int source, int target, vector<pair<int, int>>& intervals);
	bool divide_conquer_articulation_point(int source, int target, aster_result& res);
	int  divide_conquer_articulation_find(int source, int target);
	bool divide_conquer_combine(aster_result& r1,  aster_result& r2, int pivot, aster_result& comb) const;
	bool resolve_intersection_edge(int source, int target, aster_result& res);
	bool resolve_trivial_paths(int source, int target, aster_result& res);
	
	int  greedy(int source, int target);
	bool resolve_trivial_intersection(int source, int target, aster_result& res);
	int  replace_closed_nodes_w_one_edge(int source, int target, double w);

	int event_size_penalty(int eventSize) const;
	int path_distance(const path& p1, const path& p2) const;
	int edge_path_to_vertex_path(const VE& edgePath, VI& vertexPath) const;
	int find_longest_path(const aster_result& res) const;
	int find_shortest_path(const aster_result& res) const;
	int get_transcripts();
	
	string tp2v_to_string() const;
	int make_stats();
};

class aster_error : public runtime_error {
public:
    aster_error(const char* message) : runtime_error(message) {}
};

#endif