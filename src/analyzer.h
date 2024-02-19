/*
Part of Aster Assembler
(c) 2024 by Xiaofei Carl Zang, Mingfu Shao, and The Pennsylvania State University.
See LICENSE for licensing.
*/

#ifndef __ANALYZER_H__
#define __ANALYZER_H__

#include <string>
#include "splice_graph.h"
#include "hyper_set.h"
#include "equation.h"
#include "router.h"
#include "path.h"
#include "config.h"
#include "genome.h"
#include "gtf.h"
#include "aster.h"

typedef map< edge_descriptor, vector<int> > MEV;
typedef pair< edge_descriptor, vector<int> > PEV;
typedef pair< vector<int>, vector<int> > PVV;
typedef pair<PEE, int> PPEEI;
typedef map<PEE, int> MPEEI;
typedef pair<int, int> PI;
typedef map<int, int> MI;
typedef vector<int> VI;

class analyzer
{
public:
    analyzer(string);
    int write();
    int print();

private:
    genome gn;
    vector<transcript> trsts;
    vector<transcript> non_full_trsts;
    int analyze();

};

#endif