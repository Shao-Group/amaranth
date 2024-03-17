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
#include "aster_aux.hpp"

class aster
{
public:
	aster(const splice_graph &gr, const hyper_set &hs);
	int assemble();
	static int print_stats();

private:		   // static stats
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
	inline static int dnc_counter_resolve_trivial_node = 0;
	inline static int dnc_counter_resolve_intersection_edge = 0;
	int stepCount = 0;

private:
	const splice_graph& origr;			// original splice graph
	splice_graph gr;					// splice graph with modification
	hyper_set hs;						// hyper edges
	VI tp2v;							// DFS-based topologically sorted index to vertex index.
	// VI v2tp;							// vertex index to DFS-based topologically sorted index.
	map<edge_descriptor, aster_result> edgeres;
	aster_mode mode;

public:
	vector<path> paths;							// predicted paths, original v index, inclusive
    vector<transcript> trsts;					// predicted transcripts, original v index
	vector<transcript> non_full_trsts;			// predicted non full length transcripts
	bool successStatus = false;					// whether splice graph is decomposed successfully

private: 
	int init_edgeres();
	int topological_sort_vertices();
	int topological_sort_vertices_visit(int i, vector<bool>& visited, vector<int>& sorted);
	//  int topological_sort_index_edges();
	//	bool aggressive_purge_intersecting_edges();
	int prepare_graph();
	int balance_vertex(int);
	int remove_small_junctions();

	int divide_conquer();
	int divide_conquer(aster_index ai);
	bool resolve_trivial_node(aster_index ai);
	bool divide_conquer_single_vertex(aster_index ai);
	bool divide_conquer_unitig(aster_index ai);
	bool divide_conquer_abutting(aster_index ai);
	bool divide_conquer_cut_termini(aster_index ai);
	int  divide_conquer_cut_termini_find(aster_index ai, set<aster_index>& aiSubIntervals);
	bool divide_conquer_articulation_point(aster_index ai);
	int  divide_conquer_articulation_find(aster_index ai, aster_index left, aster_index right);
	// bool resolve_trivial_intersection(aster_index ai);
	bool resolve_intersection_edge(aster_index ai);
	bool greedy(aster_index ai);

	int edges_combine_consecutive_and_replace(PEEI inEdges, PEEI outEdges);
	int edge_combine_consecutive_pop_res(edge_descriptor in, edge_descriptor out, const map<pair<int, int>, vector<const path*> >& st2Path);
	bool res_combine_consecutive(aster_result& r1,  aster_result& r2, aster_result& comb) const;
	bool res_combine_parallel   (aster_result& r1,  aster_result& r2, aster_result& comb, bool frontSame = true, bool backSame = true) const;
	edge_descriptor replace_aster_index_to_one_edge(aster_index ai, double w, aster_result res);

	int non_isolated_vertex_index(aster_index& ai) const;
	bool valid_paths(aster_result res) const;
	bool valid_paths(vector<path> paths) const;
	int event_size_penalty(int eventSize) const;
	int path_distance(const path& p1, const path& p2) const;
	int edge_path_to_vertex_path(const VE& edgePath, VI& vertexPath) const;
	int find_longest_path(const aster_result& res) const;
	int find_shortest_path(const aster_result& res) const;
	int get_transcripts();
	
	string tp2v_to_string() const;
	int make_stats();
};

#endif