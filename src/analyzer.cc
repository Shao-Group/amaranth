/*
Part of Amaranth Transcript Assembler
(c) 2025 by Xiaofei Carl Zang, Mingfu Shao, and The Pennsylvania State University.
See LICENSE for licensing.
*/

#include <algorithm>
#include <cmath>
#include "util.h"
#include "amaranth.h"
#include "basic_algo.h"
#include "analyzer.h"

/* read a gtf file, make splice graphs, analyze each of them
*/
analyzer::analyzer(string _ref_file)
{
    gn = genome(_ref_file);
    analyze();
    write();
}

int analyzer::analyze()
{
    if(gn.genes.size() <= 0)
    {
        throw runtime_error("error: gtf file is empty or corrupted!");
    }

    trsts.clear();
    non_full_trsts.clear();
    
    for(int i = 0; i < gn.genes.size(); i++)
    {
        gene& gene0 = gn.genes[i];
        if(gene0.transcripts.size() <= 0) continue;
        
        gtf gtf0(gene0);
        splice_graph gr;
        gr.gid = gene0.transcripts[0].gene_id;

        gtf0.build_splice_graph(gr);
        amaranth amaranthInstance(gr, {}, true);

        trsts.insert(trsts.end(), amaranthInstance.trsts.begin(), amaranthInstance.trsts.end());
        non_full_trsts.insert(non_full_trsts.end(), amaranthInstance.non_full_trsts.begin(), amaranthInstance.non_full_trsts.end());
    }

    return 0;
}

int analyzer::print()
{
    amaranth::print_stats();
    return 0;
}

int analyzer::write()
{
    ofstream fout(output_file.c_str());
	if(fout.fail()) return 0;

	for(int i = 0; i < trsts.size(); i++)
	{
		transcript &t = trsts[i];
		t.write(fout);
	}
	fout.close();

    ofstream fout1(output_file1.c_str());
    if(fout1.fail()) return 0;
    for(int i = 0; i < non_full_trsts.size(); i++)
    {
            transcript &t = non_full_trsts[i];
            t.write(fout1);
    }
    fout1.close();

	return 0;
}