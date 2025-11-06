/*
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
Part of Scallop2
(c) 2021 by  Qimin Zhang, Mingfu Shao, and The Pennsylvania State University.
See LICENSE for licensing.
*/

#include <cstdio>
#include <cassert>
#include <sstream>
#include <map>

#include "config.h"
#include "gtf.h"
#include "genome.h"
#include "assembler.h"
#include "scallop.h"
#include "bundle.h"
#include "amaranth.h"
#include "sgraph_compare.h"
#include "super_graph.h"
#include "filter.h"

assembler::assembler()
{
    sfn = sam_open(input_file.c_str(), "r");
    hdr = sam_hdr_read(sfn);
    b1t = bam_init1();
	hid = 0;
	index = 0;
	terminate = false;
	qlen = 0;
	qcnt = 0;
}

assembler::~assembler()
{
    bam_destroy1(b1t);
    bam_hdr_destroy(hdr);
    sam_close(sfn);
}

int assembler::assemble()
{
    while(sam_read1(sfn, hdr, b1t) >= 0)
	{
		if(terminate == true) return 0;

		bam1_core_t &p = b1t->core;

		if(p.tid < 0) continue;
		if((p.flag & 0x4) >= 1) continue;										// read is not mapped
		if(remove_dup && ((p.flag & 0x400) >= 1)) continue;						// read is PCR or optical duplicate
		if((p.flag & 0x100) >= 1 && use_second_alignment == false) continue;	// secondary alignment
		if(p.n_cigar > max_num_cigar) continue;									// ignore hits with more than max-num-cigar types
		if(p.qual < min_mapping_quality) continue;							    // ignore hits with small quality
		if(p.n_cigar < 1) continue;												// should never happen

		hit ht(b1t, hid++);
		ht.set_tags(b1t);
		ht.set_strand();
		//ht.print();

		//if(ht.nh >= 2 && p.qual < min_mapping_quality) continue;
		//if(ht.nm > max_edit_distance) continue;

		//if(p.tid > 1) break;

		qlen += ht.qlen;
		qcnt += 1;

		// truncate
		if(ht.tid != bb1.tid || ht.pos > bb1.rpos + min_bundle_gap)
		{
			pool.push_back(bb1);
			bb1.clear();
		}
		if(ht.tid != bb2.tid || ht.pos > bb2.rpos + min_bundle_gap)
		{
			pool.push_back(bb2);
			bb2.clear();
		}

		// process
		process(batch_bundle_size);

		//printf("read strand = %c, xs = %c, ts = %c\n", ht.strand, ht.xs, ht.ts);

		// add hit
		if(uniquely_mapped_only == true && ht.nh != 1) continue;
		if(library_type != UNSTRANDED && ht.strand == '+' && ht.xs == '-') continue;
		if(library_type != UNSTRANDED && ht.strand == '-' && ht.xs == '+') continue;
		if(library_type != UNSTRANDED && ht.strand == '.' && ht.xs != '.') ht.strand = ht.xs;
		if(library_type != UNSTRANDED && ht.strand == '+') bb1.add_hit(ht);
		if(library_type != UNSTRANDED && ht.strand == '-') bb2.add_hit(ht);

		if (ht.umi == "" && tech == seq::SC )
		{
			if(library_type != UNSTRANDED && ht.xs == '.' && ht.spos.size() <= 0) 
			{
				ht.strand = '+';
				bb1.add_hit(ht);
			}
			if(library_type != UNSTRANDED && ht.xs == '.' && ht.spos.size() <= 0) 
			{
				ht.strand = '-';
				bb2.add_hit(ht);
			}
		}
		

		if(library_type == UNSTRANDED && ht.xs == '.' && ht.spos.size() <= 0) bb1.add_hit(ht);
		if(library_type == UNSTRANDED && ht.xs == '.' && ht.spos.size() <= 0) bb2.add_hit(ht);
		if(library_type == UNSTRANDED && ht.xs == '+') bb1.add_hit(ht);
		if(library_type == UNSTRANDED && ht.xs == '-') bb2.add_hit(ht);
	}

	pool.push_back(bb1);
	pool.push_back(bb2);
	process(0);

	amaranth::print_stats();

	assign_RPKM();

	filter ft(trsts);
	ft.merge_single_exon_transcripts();
	trsts = ft.trs;

	filter ft1(non_full_trsts);
	ft1.merge_single_exon_transcripts();
	non_full_trsts = ft1.trs;

	write();
	
	if (meta_cell_assembly) write_individual_cell();

	return 0;
}

int assembler::process(int n)
{
	if(pool.size() < n) return 0;

	for(int i = 0; i < pool.size(); i++)
	{
		bundle_base &bb = pool[i];
		bb.rm_duplicated_reads();
		bb.build_maps();
		
		// if UMI reads are too few/ low proportion, we should skip it
		bool deficient_umi_num = (bb.umi_reads < min_umi_reads_bundle);
		bool deficient_umi_ratio = ((float(bb.umi_reads) / float(bb.hits.size())) < min_umi_ratio_bundle);
		if(both_umi_support)
		{
			if(deficient_umi_num || deficient_umi_ratio) continue;
		}
		else
		{
			if(deficient_umi_num && deficient_umi_ratio) continue;
		}

		int cnt1 = 0;
		int cnt2 = 0;
		for(int k = 0; k < bb.hits.size(); k++)
		{
			//counts += (1 + bb.hits[k].spos.size());
			if(bb.hits[k].spos.size() >= 1) cnt1 ++;
			else cnt2++;
		}

		if(cnt1 + cnt2 < min_num_hits_in_bundle) continue;
		//if(cnt1 < 5 && cnt1 * 2 + cnt2 < min_num_hits_in_bundle) continue;
		if(bb.tid < 0) continue;

		char buf[1024];
		strcpy(buf, hdr->target_name[bb.tid]);
		bb.chrm = string(buf);

		transcript_set ts1(bb.chrm, 0.9);		// full-length set
		transcript_set ts2(bb.chrm, 0.9);		// non-full-length set

		bundle bd(bb);

		bd.build(1, true);
		bd.print(index++);
		assemble(bd.gr, bd.hs, ts1, ts2);

		bd.build(2, true);
		bd.print(index++);
		assemble(bd.gr, bd.hs, ts1, ts2);

		int sdup = assemble_duplicates / 1 + 1;
		int mdup = assemble_duplicates / 2 + 0;

		vector<transcript> gv1 = ts1.get_transcripts(sdup, mdup);
		vector<transcript> gv2 = ts2.get_transcripts(sdup, mdup);

		for(int k = 0; k < gv1.size(); k++)
		{
			if(gv1[k].exons.size() >= 2) gv1[k].coverage /= (1.0 * assemble_duplicates);
		}
		for(int k = 0; k < gv2.size(); k++) 
		{
			if(gv2[k].exons.size() >= 2) gv2[k].coverage /= (1.0 * assemble_duplicates);
		}
		filter ft1(gv1);
		if (use_filter)
		{
			ft1.filter_length_coverage();
			ft1.remove_nested_transcripts();
		}
		if(ft1.trs.size() >= 1) trsts.insert(trsts.end(), ft1.trs.begin(), ft1.trs.end());

		filter ft2(gv2);
		if (use_filter)
		{
			ft2.filter_length_coverage();
			ft2.remove_nested_transcripts();
		}
		if(ft2.trs.size() >= 1) non_full_trsts.insert(non_full_trsts.end(), ft2.trs.begin(), ft2.trs.end());
	}
	pool.clear();
	return 0;
}

int assembler::assemble(const splice_graph &gr0, const hyper_set &hs0, transcript_set &ts1, transcript_set &ts2)
{
	super_graph sg(gr0, hs0);
	sg.build();

	/*
	vector<transcript> gv;
	vector<transcript> gv1;
	*/

	for(int k = 0; k < sg.subs.size(); k++)
	{
		splice_graph &gr = sg.subs[k];
		hyper_set &hs = sg.hss[k];

		if(determine_regional_graph(gr) == true) continue;
		if(gr.num_edges() <= 0) continue;

		for(int r = 0; r < assemble_duplicates; r++)
		{
			string gid = "gene." + tostring(index) + "." + tostring(k) + "." + tostring(r);
			gr.gid = gid;
			
			if (algo == "amaranth" && gr.num_vertices() < 1000)
			{
				amaranth amaranthInstance(gr, hs, true);
				if(verbose >= 2)
				{
					printf("assembly with r = %d, total %lu transcripts, run 1:\n", r, amaranthInstance.trsts.size());
					for(int i = 0; i < amaranthInstance.trsts.size(); i++) amaranthInstance.trsts[i].write(cout);
				}
	
				for(int i = 0; i < amaranthInstance.trsts.size(); i++)
				{
					ts1.add(amaranthInstance.trsts[i], 1, 0, TRANSCRIPT_COUNT_ADD_COVERAGE_MIN, TRANSCRIPT_COUNT_ADD_COVERAGE_ADD);
				}
				for(int i = 0; i < amaranthInstance.non_full_trsts.size(); i++)
				{
					ts2.add(amaranthInstance.non_full_trsts[i], 1, 0, TRANSCRIPT_COUNT_ADD_COVERAGE_MIN, TRANSCRIPT_COUNT_ADD_COVERAGE_ADD);
				}
			}
			else
			{
				scallop scallopInstance(gr, hs, true);
				if(verbose >= 2)
				{
					printf("assembly with r = %d, total %lu transcripts, run 1:\n", r, scallopInstance.trsts.size());
					for(int i = 0; i < scallopInstance.trsts.size(); i++) scallopInstance.trsts[i].write(cout);
				}
				for(int i = 0; i < scallopInstance.trsts.size(); i++)
				{
					ts1.add(scallopInstance.trsts[i], 1, 0, TRANSCRIPT_COUNT_ADD_COVERAGE_MIN, TRANSCRIPT_COUNT_ADD_COVERAGE_ADD);
				}
				for(int i = 0; i < scallopInstance.non_full_trsts.size(); i++)
				{
					ts2.add(scallopInstance.non_full_trsts[i], 1, 0, TRANSCRIPT_COUNT_ADD_COVERAGE_MIN, TRANSCRIPT_COUNT_ADD_COVERAGE_ADD);
				}
			}
		}
	}

	return 0;
}

bool assembler::determine_regional_graph(splice_graph &gr)
{
	bool all_regional = true;
	for(int i = 1; i < gr.num_vertices() - 1; i++)
	{
		if(gr.get_vertex_info(i).regional == false) all_regional = false;
		if(all_regional == false) break;
	}
	return all_regional;
}

int assembler::assign_RPKM()
{
	double factor = 1e9 / qlen;
	for(int i = 0; i < trsts.size(); i++)
	{
		trsts[i].assign_RPKM(factor);
	}
	return 0;
}

int assembler::write()
{
	ofstream fout(output_file.c_str());
	if(fout.fail()) return 0;
	for(int i = 0; i < trsts.size(); i++)
	{
		transcript &t = trsts[i];
		t.write(fout);
	}
	fout.close();

	if (output_file1 != "")
	{
		ofstream fout1(output_file1.c_str());
		if(!fout1.fail())
		{
			for(int i = 0; i < non_full_trsts.size(); i++)
			{
				transcript &t = non_full_trsts[i];
				t.write(fout1);
			}
			fout1.close();
		}
	}

	// Write features if output_feat is specified
	if(output_feat != "")
	{
		write_features();
	}

	return 0;
}

int assembler::write_features()
{
	if(output_feat == "") return 0;

	ofstream fout(output_feat.c_str());
	if(fout.fail())
	{
		cerr << "Error: Cannot open feature output file " << output_feat << endl;
		return 0;
	}

	// Combine all transcripts for feature output
	vector<transcript*> all_trsts;
	for(int i = 0; i < trsts.size(); i++)
	{
		all_trsts.push_back(&trsts[i]);
	}
	for(int i = 0; i < non_full_trsts.size(); i++)
	{
		all_trsts.push_back(&non_full_trsts[i]);
	}

	if(all_trsts.empty()) return 0;

	// Get all feature keys from the first transcript
	vector<string> feature_keys;
	for(const auto &kv : all_trsts[0]->features)
	{
		feature_keys.push_back(kv.first);
	}
	sort(feature_keys.begin(), feature_keys.end());  // Sort for consistent column order

	// Write header: basic info first, then features
	fout << "chr,start,end,strand,transcript_id,gene_id";
	for(const auto &key : feature_keys)
	{
		fout << "," << key;
	}
	fout << endl;

	// Write features for each transcript
	for(int i = 0; i < all_trsts.size(); i++)
	{
		transcript *t = all_trsts[i];

		// Assert that all transcripts have the same keys
		assert(t->features.size() == all_trsts[0]->features.size());
		for(const auto &key : feature_keys)
		{
			assert(t->features.find(key) != t->features.end());
		}

		// Get transcript bounds
		PI32 bounds = t->get_bounds();

		// Write basic transcript information
		fout << t->seqname << ",";           // chr
		fout << bounds.first + 1 << ",";     // start (1-based)
		fout << bounds.second << ",";        // end
		fout << t->strand << ",";            // strand
		fout << t->transcript_id << ",";     // transcript_id
		fout << t->gene_id;                  // gene_id

		// Write feature values
		for(const auto &key : feature_keys)
		{
			fout << ",";
			// Use std::visit to handle variant types
			std::visit([&fout](const auto& value) {
				fout << value;
			}, t->features[key]);
		}
		fout << endl;
	}

	fout.close();
	return 0;
}

int assembler::compare(splice_graph &gr, const string &file, const string &texfile)
{
	if(file == "") return 0;

	genome g(file);
	if(g.genes.size() <= 0) return 0;

	gtf gg(g.genes[0]);

	splice_graph gt;
	gg.build_splice_graph(gt);

	sgraph_compare sgc(gt, gr);
	sgc.compare(texfile);

	return 0;
}
