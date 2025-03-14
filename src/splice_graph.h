/*
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
Part of Scallop2
(c) 2021 by  Qimin Zhang, Mingfu Shao, and The Pennsylvania State University.
See LICENSE for licensing.
*/

#ifndef __SPLICE_GRAPH_H__
#define __SPLICE_GRAPH_H__

#include "directed_graph.h"
#include "vertex_info.h"
#include "edge_info.h"
#include "path.h"
#include "gene.h"

#include <map>
#include <cassert>

#define SMIN 0.00001

using namespace std;

typedef map<edge_descriptor, edge_info> MEIF;
typedef pair<edge_descriptor, edge_info> PEIF;

class splice_graph : public directed_graph
{
public:
	splice_graph();
	splice_graph(const splice_graph &gr);
	virtual ~splice_graph();

public:
	string chrm;
	string gid;
	char strand;

	vector<double> vwrt;
	vector<vertex_info> vinf;
	MED ewrt;
	MEIF einf;

public:
	// get and set properties
	double get_vertex_weight(int v) const;
	double get_edge_weight(edge_base *e) const;
	vertex_info get_vertex_info(int v) const;
	edge_info get_edge_info(edge_base *e) const;

	int set_vertex_weight(int v, double w);
	int set_vertex_info(int v, const vertex_info &vi);
	int set_edge_weight(edge_base *e, double w);
	int set_edge_info(edge_base *e, const edge_info &ei);

	MED get_edge_weights() const;
	vector<double> get_vertex_weights() const;
	int set_edge_weights(const MED & med);
	int set_vertex_weights(const vector<double> &v);

	edge_descriptor max_out_edge(int v);
	edge_descriptor max_in_edge(int v);
	double get_in_weights(int v);
	double get_out_weights(int v);
	double get_max_in_weight(int v);
	double get_max_out_weight(int v);

	int count_junctions();

	// modify the splice_graph
	int clear();
	int copy(const splice_graph &gr, MEE &x2y, MEE &y2x);
	bool refine_splice_graph();

	// read, write, and simulate splice graph
	int build(const string &file);
	int write(const string &file) const;
	int simulate(int nv, int ne, int mf);

	// analysis the structure of splice graph
	long compute_num_paths();
	long compute_num_paths(int a, int b);
	long compute_num_paths(int a, int b, int max);
	int compute_decomp_paths();
	bool check_fully_connected();
	int compute_independent_subgraphs();
	double compute_average_vertex_weight();
	double compute_average_edge_weight();

	// algorithms with weight contraints
	edge_descriptor compute_maximum_edge_w();
	int bfs_w(int s, double w, vector<int> &v, VE &b);
	int compute_shortest_path_w(int s, int t, double w);
	int compute_shortest_path_w(int s, int t, double w, VE &p);
	int compute_closest_path(int s, vector<double> &d);
	int compute_closest_path(int s, vector<double> &d, vector<int> &b);
	int compute_closest_path_reverse(int t, vector<double> &d);
	int compute_closest_path_reverse(int t, vector<double> &d, vector<int> &b);
	double compute_maximum_path_w(VE &p);
	double compute_maximum_st_path_w(VE &p, int s, int t);
	double compute_minimum_weight(const VE &p);

	// determine optimal path
	bool compute_optimal_path(VE &p);

	// rounding all weights to integers
	int round_weights();
	int locate(int v);
	int locate_vertex(int32_t p);
	int locate_vertex(int32_t p, int a, int b);

	// draw and print
	int draw(const string &file, string label = "");
	int draw(const string &file, const MIS &mis, const MES &mes, double len, string label = "");
	int draw(const string &file, const MIS &mis, const MES &mes, double len, const vector<int> &tp, string label = "");
	int graphviz(const string &file, string label = "");
	int graphviz(const string &file, const MIS &mis, const MES &mes, double len, string label = "");
	int graphviz(const string &file, const MIS &mis, const MES &mes, double len, const vector<int> &tp, string label = "");
	int print_nontrivial_vertices();
	int print_weights();
	int print();

	// output transcripts
	int output_transcripts(ofstream &fout, const vector<path> &p) const;
	int output_transcripts(vector<transcript> &v, const vector<path> &p) const;
	int output_transcripts1(vector<transcript> &v, vector<transcript> &v1, const vector<path> &p) const;
	int output_transcript(ofstream &fout, const path &p, const string &tid) const;
	int output_transcript(transcript &trst, const path &p, const string &tid) const;
};

#endif
