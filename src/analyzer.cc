/*
Part of Aster Transcript Assembler
(c) 2024 by Xiaofei Carl Zang, Mingfu Shao, and The Pennsylvania State University.
See LICENSE for licensing.
*/

#include <algorithm>
#include <cmath>
#include "util.h"
#include "aster.h"
#include "basic_algo.h"
#include "analyzer.h"

/* read a gtf file, make splice graphs, analyze each of them
*/
analyzer::analyzer(string _ref_file)
{
    genome g(_ref_file);
    analyze();
}

int analyzer::analyze()
{
    if(gn.genes.size() <= 0)
    {
        throw runtime_error("error: gtf file is empty or corrupted!");
    }
    for(int i = 0; i < gn.genes.size(); i++)
    {
        gtf gg(gn.genes[i]);
        splice_graph gr;
        gg.build_splice_graph(gr);
        aster asterInstance(gr, {});
    }

    return 0;
}

int analyzer::print()
{
    aster::print_stats();
    return 0;
}