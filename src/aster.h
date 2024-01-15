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


class aster
{
public:
	aster(const splice_graph &gr, const hyper_set &hs, bool random_ordering = false);
	int assemble();
	static int print_stats();

// static stats
private:
	inline static int num_graph = 0;
	inline static int num_intron = 0;
	inline static int num_exon = 0;
	inline static int num_overlapping_intron_count = 0;
	inline static int num_overlapping_intron_pair = 0;

public:
	splice_graph gr;					// splice graph
	hyper_set hs;						// hyper edges
	MEI e2i;							// edge map, from edge to index, sorted
	VE i2e;								// edge map, from index to edge, sorted
    vector<transcript> trsts;			// predicted transcripts
	vector<transcript> non_full_trsts;		// predicted non full length transcripts

private:
	bool random_ordering;				// whether using random ordering
	int round;							// iteration
	vector<path> paths;					// predicted paths

private: 
	int topological_sort_vertices();
	int topological_sort_index_edges();
	int make_stats();
};

class astron
{
public:
	astron(const aster*, const vector<int>& _canons, const vector<int>& _illegal, const vector<int>& _alts = {});

private:
	const aster* as;
    vector<int> canons;		// canonical events
	vector<int> illegals;		// illegal events
    vector<int> alternatives;	// alternative events

public:
	int dist;
	vector<path> paths;

private:
	int classify();
	int divide_and_conquer();
	int dnc_combine(const vector<path> subpaths, int eventOfConcern);
	int collect_trivial_path();

	// help functions	
	int event_size_penalty(int eventSize);
};

#endif