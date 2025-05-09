/*
Part of Scallop Transcript Assembler
(c) 2017 by Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
Part of Coral
(c) 2019 by Mingfu Shao, The Pennsylvania State University.
Part of Scallop2
(c) 2021 by  Qimin Zhang, Mingfu Shao, and The Pennsylvania State University.
See LICENSE for licensing.
*/

#include "config.h"
#include <cstdlib>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <cstring>

using namespace std;

//// parameters
// for bam file and reads
int min_flank_length = 3;
int max_num_cigar = 1000;
int max_edit_distance = 10;
int32_t min_bundle_gap = 100;		
int min_num_hits_in_bundle = 5;	
int min_num_splices_in_bundle = 15;	// not used; accept bundle if #hits with splices is at least this number
uint32_t min_mapping_quality = 1;
int32_t min_splice_boundary_hits = 1;
bool use_second_alignment = false;
bool uniquely_mapped_only = false;
int library_type = EMPTY;

// for preview
int max_preview_reads = 2000000;
int max_preview_spliced_reads = 50000;
int min_preview_spliced_reads = 10000;
int max_preview_umi_reads = 10000;
double preview_infer_ratio = 0.85;
bool preview_only = false;
double insertsize_ave = 300;
double insertsize_std = 50;
int insertsize_median = -1;
int insertsize_low = -1;
int insertsize_high = -1;
double insertsize_low_percentile = 0.005;
double insertsize_high_percentile = 0.998;
vector<double> insertsize_profile;

// for bridging
double min_bridging_score = 0.5;
int max_num_path_nodes = 10000;
int dp_solution_size = 10;
int dp_stack_size = 5;
bool use_overlap_scoring = false;
int32_t max_clustering_flank = 30;
int32_t flank_tiny_length = 10;
double flank_tiny_ratio = 0.4;

// for identifying subgraphs
int32_t min_subregion_gap = 3;
int32_t min_subregion_len = 15;
int32_t min_subregion_max = 3;
double min_subregion_ave = 1.5;

// for amaranth
amaranth_mode amaranthMode = amaranth_mode::ASSEMBLER;
amaranth_strategy amaranthStrategy = amaranth_strategy::HEAVY;
seq tech = seq::UNKNOWN;
bool doesFastDnC = false;

// for revising/decomposing splice graph
double min_guaranteed_edge_weight = 0.01;
double min_surviving_edge_weight = 1.5;
double max_intron_contamination_coverage = 2.0;
double max_decompose_error_ratio[7] = {0.33, 0.05, 0.0, 0.25, 0.30, 0.0, 1.1};

// for selecting paths
double min_transcript_coverage = 1.5;
double min_transcript_coverage_ratio = 0.005;
double min_single_exon_coverage = 20;
double min_transcript_numreads = 10;
int min_transcript_length_base = 150;
int min_transcript_length_increase = 50;
int min_exon_length = 20;
int max_num_exons = 1000;

// for subsetsum and router
int max_dp_table_size = 10000;
int min_router_count = 1;

// for simulation
int simulation_num_vertices = 0;
int simulation_num_edges = 0;
int simulation_max_edge_weight = 0;

// input and output
string algo = "amaranth";
string input_file;
string ref_file;
string output_file;
string output_file1;

// umi & hybrid parameters
int min_umi_reads_bundle = 1;
double min_umi_ratio_bundle = 0;
bool both_umi_support = false;
int min_umi_reads_start_exon = 1;

// filtering & retention
bool remove_retained_intron = true;
int remove_dup = 0;              		// 0: nothing, 1: by algin+cigar, 2: algin+cigar.w.S
bool use_filter = true;

// for controling
bool output_tex_files = false;
bool output_graphviz_files = false;
string fixed_gene_name = "";
int batch_bundle_size = 100;
int verbose = 1;
int assemble_duplicates = 1; //TODO:
string version = "v0.1";

int parse_arguments(int argc, const char ** argv)
{
	for(int i = 1; i < argc; i++)
	{
		// necessary ones
		if(string(argv[i]) == "-i")
		{
			input_file = string(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "-o")
		{
			output_file = string(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "-f")
		{
			output_file1 = string(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--transcript_fragments")
		{
			output_file1 = string(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "-r")
		{
			ref_file = string(argv[i + 1]);
			i++;
		}

		// internal use
		else if(string(argv[i]) == "-g")
		{
			fixed_gene_name = string(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "-t")
		{
			output_tex_files = true;
		}
		else if(string(argv[i]) == "-z")
		{
			output_graphviz_files = true;
		}

		// user specified
		else if(string(argv[i]) == "--version")
		{
			printf("%s\n", version.c_str());
			exit(0);
		}
		else if(string(argv[i]) == "--help")
		{
			print_copyright();
			print_help();
			printf("\n");
			exit(0);
		}
		else if(string(argv[i]) == "--min_flank_length")
		{
			min_flank_length = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--max_num_cigar")
		{
			max_num_cigar = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--max_edit_distance")
		{
			max_edit_distance = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--min_bundle_gap")
		{
			min_bundle_gap = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--min_num_hits_in_bundle")
		{
			min_num_hits_in_bundle = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--min_num_splices_in_bundle")
		{
			min_num_splices_in_bundle = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--min_mapping_quality")
		{
			min_mapping_quality = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--min_splice_boundary_hits")
		{
			min_splice_boundary_hits = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--max_preview_spliced_reads")
		{
			max_preview_spliced_reads = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--min_preview_spliced_reads")
		{
			min_preview_spliced_reads = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--preview")
		{
			preview_only = true;
		}
		else if(string(argv[i]) == "--max_preview_reads")
		{
			max_preview_reads = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--preview_infer_ratio")
		{
			preview_infer_ratio = atof(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--min_subregion_gap")
		{
			min_subregion_gap = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--min_subregion_len")
		{
			min_subregion_len = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--min_subregion_ave")
		{
			min_subregion_ave = atof(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--min_subregion_max")
		{
			min_subregion_max = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--min_surviving_edge_weight")
		{
			min_surviving_edge_weight = atof(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--max_intron_contamination_coverage")
		{
			max_intron_contamination_coverage = atof(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--min_transcript_coverage")
		{
			min_transcript_coverage = atof(argv[i + 1]);
			i++;
			if(fabs(min_transcript_coverage - 1.0) < 0.01) min_transcript_coverage = 1.01;
		}
		else if(string(argv[i]) == "--min_transcript_coverage_ratio")
		{
			min_transcript_coverage_ratio = atof(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--min_single_exon_coverage")
		{
			min_single_exon_coverage = atof(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--min_transcript_numreads")
		{
			min_transcript_numreads = atof(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--min_transcript_length_base")
		{
			min_transcript_length_base = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--min_transcript_length_increase")
		{
			min_transcript_length_increase = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--min_exon_length")
		{
			min_exon_length = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--max_num_exons")
		{
			max_num_exons = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--max_dp_table_size")
		{
			max_dp_table_size = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--min_router_count")
		{
			min_router_count = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--amaranthMode" || string(argv[i]) == "-m")
		{
			string s(argv[i + 1]);
			if(s == "REF") amaranthMode = amaranth_mode::REF;
			else if(s == "STAT") amaranthMode = amaranth_mode::STAT_ONLY;
			else if(s == "MINI") amaranthMode = amaranth_mode::MINI;
			else if(s == "ASSEMBLY") amaranthMode = amaranth_mode::ASSEMBLER;
			else throw runtime_error("received unknown --amaranthMode value: " + s + " (must be REF/STAT/MINI/ASSEMBLY)");
			i++;
		}
		else if(string(argv[i]) == "--amaranthStrategy" || string(argv[i]) == "-s")
		{
			string s(argv[i + 1]);
			if(s == "LONG") amaranthStrategy = amaranth_strategy::LONG;
			else if(s == "SHORT") amaranthStrategy = amaranth_strategy::SHORT;
			else if(s == "HEAVY") amaranthStrategy = amaranth_strategy::HEAVY;
			else throw runtime_error("received unknown --amaranthStrategy value: " + s + " (must be LONG/SHORT/HEAVY)");
			i++;
		}
		else if(string(argv[i]) == "--tech")
		{
			string s(argv[i + 1]);
			if(s == "SC") tech = seq::SC;
			else if(s == "BULK") tech = seq::BULK;
			else if(s == "UNKNOWN") tech = seq::UNKNOWN; 
			else throw runtime_error("received unknown --tech value: " + s + " (must be SC/BULK)");
			i++;
		}
		else if(string(argv[i]) == "--fast")
		{
			doesFastDnC = true;
		}
		else if(string(argv[i]) == "--no-filter")
		{
			use_filter = false;
		}
		else if(string(argv[i]) == "--use-filter")
		{
			use_filter = true;
		}
		else if(string(argv[i]) == "--remove-pcr-duplicates")
		{
			remove_dup = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--remove-reteind-intron")
		{
			remove_retained_intron = true;
		}
		else if(string(argv[i]) == "--no-remove-reteind-intron")
		{
			remove_retained_intron = false;
		}
		else if(string(argv[i]) == "--min_umi_reads_bundle")
		{
			min_umi_reads_bundle = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--min_umi_ratio_bundle")
		{
			min_umi_ratio_bundle = atof(argv[i + 1]);
			if (min_umi_ratio_bundle < 0 || min_umi_ratio_bundle > 1)
			{
				throw runtime_error("min_umi_ratio_bundle must be between 0 and 1.");
			}
			i++;
		}
		else if(string(argv[i]) == "--both_umi_support")
		{
			string s(argv[i + 1]);
			if(s == "true" || s == "True" || s == "TRUE" || s == "T") both_umi_support = true;
			else if(s == "false" || s == "False" || s == "FALSE" || s == "F") both_umi_support = false;
			else throw runtime_error("received unknown --both_umi_support value: " + s + " (must be true/false)");
			i++;
		}
		else if(string(argv[i]) == "--min_umi_reads_start_exon")
		{
			min_umi_reads_start_exon = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--library_type")
		{
			string s(argv[i + 1]);
			if(s == "empty") library_type = EMPTY;
			else if(s == "unstranded") library_type = UNSTRANDED;
			else if(s == "first") library_type = FR_FIRST;
			else if(s == "second") library_type = FR_SECOND;
			else throw runtime_error("received unknown --library_type value: " + s + " (must be empty/unstranded/first/second)");
			i++;
		}
		else if(string(argv[i]) == "--use_second_alignment")
		{
			string s(argv[i + 1]);
			if(s == "true") use_second_alignment = true;
			else use_second_alignment = false;
			i++;
		}
		else if(string(argv[i]) == "--uniquely_mapped_only")
		{
			string s(argv[i + 1]);
			if(s == "true") uniquely_mapped_only = true;
			else uniquely_mapped_only = false;
			i++;
		}
		else if(string(argv[i]) == "--verbose")
		{
			verbose = atoi(argv[i + 1]);
			i++;
		}
		// else if(string(argv[i]) == "--assemble_duplicates")
		// {
		// 	assemble_duplicates = atoi(argv[i + 1]);
		// 	i++;
		// }
		else if(string(argv[i]) == "--batch_bundle_size")
		{
			batch_bundle_size = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--min_bridging_score")
		{
			min_bridging_score = atof(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--dp_solution_size")
		{
			dp_solution_size = atof(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--dp_stack_size")
		{
			dp_stack_size = atof(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--max_clustering_flank")
		{
			max_clustering_flank = atof(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--flank_tiny_length")
		{
			flank_tiny_length = atof(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--flank_tiny_ratio")
		{
			flank_tiny_ratio = atof(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--insertsize_median")
		{
			insertsize_median = atof(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--insertsize_low")
		{
			insertsize_low = atof(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--insertsize_high")
		{
			insertsize_high = atof(argv[i + 1]);
			i++;
		}
		else
		{
			throw runtime_error("Unknown argument is provided: " + string(argv[i]));
		}
	}

	if(min_surviving_edge_weight < 0.1 + min_transcript_coverage) 
	{
		min_surviving_edge_weight = 0.1 + min_transcript_coverage;
		if(min_surviving_edge_weight > 10.0) min_surviving_edge_weight = 10.0;
	}

	// verify arguments
	if(amaranthMode != amaranth_mode::REF)
	{
		if(input_file == "")  
		{
			throw runtime_error("error: input-file is missing.");
		}	
		if(output_file == "" && preview_only == false)  
		{
			throw runtime_error("error: output-file is missing.");
		}
		
	}
	if(amaranthMode == amaranth_mode::REF)
	{
		if(ref_file == "") 
		{
			throw runtime_error("error: amaranthMode is ref, but reference-file is missing.");
		}
	}


	return 0;
}

int print_parameters()
{
	printf("parameters:\n");

	// for bam file and reads
	printf("min_flank_length = %d\n", min_flank_length);
	printf("max_num_cigar = %d\n", max_num_cigar);
	printf("max_edit_distance = %d\n", max_edit_distance);
	printf("min_bundle_gap = %d\n", min_bundle_gap);
	printf("min_num_hits_in_bundle = %d\n", min_num_hits_in_bundle);
	printf("min_mapping_quality = %d\n", min_mapping_quality);
	printf("min_splice_boundary_hits = %d\n", min_splice_boundary_hits);

	// for preview
	printf("preview_only = %c\n", preview_only ? 'T' : 'F');
	printf("max_preview_reads = %d\n", max_preview_reads);
	printf("max_preview_spliced_reads = %d\n", max_preview_spliced_reads);
	printf("min_preview_spliced_reads = %d\n", min_preview_spliced_reads);
	printf("preview_infer_ratio = %.3lf\n", preview_infer_ratio);

	// for identifying subgraphs
	printf("min_subregion_gap = %d\n", min_subregion_gap);
	printf("min_subregion_len = %d\n", min_subregion_len);
	printf("min_subregion_max = %d\n", min_subregion_max);
	printf("min_subregion_ave = %.2lf\n", min_subregion_ave);

	// for splice graph
	printf("max_intron_contamination_coverage = %.2lf\n", max_intron_contamination_coverage);
	printf("min_surviving_edge_weight = %.2lf\n", min_surviving_edge_weight);
	printf("min_transcript_coverage = %.2lf\n", min_transcript_coverage);
	printf("min_transcript_coverage_ratio = %.2lf\n", min_transcript_coverage_ratio);
	printf("min_single_exon_coverage = %.2lf\n", min_single_exon_coverage);
	printf("min_transcript_numreads = %.2lf\n", min_transcript_numreads);
	printf("min_transcript_length_base = %d\n", min_transcript_length_base);
	printf("min_transcript_length_increase = %d\n", min_transcript_length_increase);
	printf("max_num_exons = %d\n", max_num_exons);

	// for subsetsum and router
	printf("max_dp_table_size = %d\n", max_dp_table_size);
	printf("min_router_count = %d\n", min_router_count);

	// for simulation
	printf("simulation_num_vertices = %d\n", simulation_num_vertices);
	printf("simulation_num_edges = %d\n", simulation_num_edges);
	printf("simulation_max_edge_weight = %d\n", simulation_max_edge_weight);

	// for input and output
	printf("algo = %s\n", algo.c_str());
	printf("input_file = %s\n", input_file.c_str());
	printf("ref_file = %s\n", ref_file.c_str());
	printf("output_file = %s\n", output_file.c_str());
	printf("output_file1 = %s\n", output_file1.c_str());

	// for controling
	printf("library_type = %d\n", library_type);
	printf("output_tex_files = %c\n", output_tex_files ? 'T' : 'F');
	printf("output_graphviz_files = %c\n", output_graphviz_files ? 'T' : 'F');
	printf("fixed_gene_name = %s\n", fixed_gene_name.c_str());
	printf("use_second_alignment = %c\n", use_second_alignment ? 'T' : 'F');
	printf("uniquely_mapped_only = %c\n", uniquely_mapped_only ? 'T' : 'F');
	printf("verbose = %d\n", verbose);
	printf("batch_bundle_size = %d\n", batch_bundle_size);

	printf("\n");

	return 0;
}

int print_command_line(int argc, const char ** argv)
{
	printf("command line: ");
	for(int i = 0; i < argc; i++)
	{
		printf("%s ", argv[i]);
	}
	printf("\n");
	return 0;
}

int print_help()
{
	printf("\n");
	printf("Usage: amaranth -i <bam-file> -o <gtf-file> [options]\n");
	printf("\n");
	printf("Options:\n");
	printf(" %-42s  %s\n", "--help",  "print usage of Scallop and exit");
	printf(" %-42s  %s\n", "--version",  "print current version of Scallop and exit");
	printf(" %-42s  %s\n", "--preview",  "determine fragment-length-range and library-type and exit");
	printf(" %-42s  %s\n", "--verbose <0, 1, 2>",  "0: quiet; 1: one line for each graph; 2: with details, default: 1");
	printf(" %-42s  %s\n", "-f/--transcript_fragments <filename>",  "file to which the assembled non-full-length transcripts will be written to");
	printf(" %-42s  %s\n", "--library_type <first, second, unstranded>",  "library type of the sample, default: unstranded");
	// printf(" %-42s  %s\n", "--assemble_duplicates <integer>",  "the number of consensus runs of the decomposition, default: 10");
	printf(" %-42s  %s\n", "--min_transcript_coverage <float>",  "minimum coverage required for a multi-exon transcript, default: 1.5");
	printf(" %-42s  %s\n", "--min_single_exon_coverage <float>",  "minimum coverage required for a single-exon transcript, default: 20");
	printf(" %-42s  %s\n", "--min_transcript_length_increase <integer>",  "default: 50");
	printf(" %-42s  %s\n", "--min_transcript_length_base <integer>",  "default: 150, minimum length of a transcript would be");
	printf(" %-42s  %s\n", "",  "--min_transcript_length_base + --min_transcript_length_increase * num-of-exons");
	printf(" %-42s  %s\n", "--min_mapping_quality <integer>",  "ignore reads with mapping quality less than this value, default: 1");
	printf(" %-42s  %s\n", "--max_num_cigar <integer>",  "ignore reads with CIGAR size larger than this value, default: 1000");
	printf(" %-42s  %s\n", "--min_bundle_gap <integer>",  "minimum distances required to start a new bundle, default: 100");
	printf(" %-42s  %s\n", "--min_num_hits_in_bundle <integer>",  "minimum number of reads required in a gene locus, default: 5");
	printf(" %-42s  %s\n", "--min_flank_length <integer>",  "minimum match length in each side for a spliced read, default: 3");

	printf(" %-42s  %s\n", "--use-filter",  "use filtering to select subpaths before final assembly, default: use-filter");
	printf(" %-42s  %s\n", "--no-filter",   "disable filtering, use all subpaths in final assembly,  default: use-filter");
	printf(" %-42s  %s\n", "--remove-pcr-duplicates <int>",     "remove PCR duplicates using strategy: 0,1, or 2, default: 0");
	printf(" %-42s  %s\n", "--no-remove-pcr-duplicates",  "not remove PCR duplicates in the input bam file, default: not-remove");

	// umi support settings
	printf(" %-42s  %s\n", "--min_umi_reads_bundle <integer>", "minimum number of UMI reads required in a bundle, default: 0");
	printf(" %-42s  %s\n", "--min_umi_ratio_bundle <float>", "minimum ratio of UMI reads required in a bundle, default: 0");
	printf(" %-42s  %s\n", "--both_umi_support <true|false>", "require satisfactory UMI support for [both/either] condition, default: false(either)");

	// amaranth - mode
	printf(" %-42s  %s\n", "-m/--amaranthMode <REF|STAT|MINI|ASSEMBLY>", "set AMARANTH operation mode (default: ASSEMBLY):");
	printf(" %-42s  %s\n", "", "    REF: reference-guided assembly mode");
	printf(" %-42s  %s\n", "", "    STAT: statistical analysis only mode");
	printf(" %-42s  %s\n", "", "    MINI: minimal assembly mode");
	printf(" %-42s  %s\n", "", "    ASSEMBLY: full assembly mode");

	// amaranth - path selection modes
	printf(" %-42s  %s\n", "-s/--amaranthStrategy <LONG|SHORT|HEAVY>", "set AMARANTH path selection strategy (default: HEAVY):");
	printf(" %-42s  %s\n", "", "    LONG: use longest subpath as anchor");
	printf(" %-42s  %s\n", "", "    SHORT: use shortest subpath as anchor");
	printf(" %-42s  %s\n", "", "    HEAVY: use highest coverage subpath as anchor");

	// amaranth - sequencing technology
	printf(" %-42s  %s\n", "--tech <SC|BULK>", "set sequencing technology (default: BULK):");
	printf(" %-42s  %s\n", "", "    SC: single-cell RNA-seq");
	printf(" %-42s  %s\n", "", "    BULK: bulk RNA-seq");

	return 0;
}

int print_copyright()
{
	printf("Amaranth assembler %s (c) 2024 Xiaofei Carl Zang, and Mingfu Shao, The Pennsylvania State University\n", version.c_str());
	return 0;
}

