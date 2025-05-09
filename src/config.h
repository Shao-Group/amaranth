/*
Part of Scallop Transcript Assembler
(c) 2017 by Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
Part of Coral
(c) 2019 by Mingfu Shao, The Pennsylvania State University.
Part of Scallop2
(c) 2021 by  Qimin Zhang, Mingfu Shao, and The Pennsylvania State University.
See LICENSE for licensing.
*/

#ifndef __CONFIG_H__
#define __CONFIG_H__

#include "util.h"
#include <stdint.h>
#include <map>
#include <sstream>

using namespace std;

//// constants
#define START_BOUNDARY 1
#define END_BOUNDARY 2
#define LEFT_SPLICE 3
#define RIGHT_SPLICE 4
#define LEFT_RIGHT_SPLICE 5
#define MIDDLE_CUT 6

#define TRIVIAL 0
#define NORMAL 1

// five types for decomposition
#define SMALLEST_EDGE 0
#define NEGLIGIBLE_EDGE 1
#define SPLITTABLE_SIMPLE 2
#define SPLITTABLE_HYPER 3
#define UNSPLITTABLE_SINGLE 4
#define UNSPLITTABLE_MULTIPLE 5
#define TRIVIAL_VERTEX 6

#define EMPTY -1
#define UNSTRANDED 0
#define FR_FIRST 1
#define FR_SECOND 2

#define EMPTY_VERTEX -9

#define TRANSCRIPT_COUNT_ADD_COVERAGE_ADD 1
#define TRANSCRIPT_COUNT_ADD_COVERAGE_NUL 2
#define TRANSCRIPT_COUNT_ADD_COVERAGE_MAX 3
#define TRANSCRIPT_COUNT_ADD_COVERAGE_MIN 4

//// parameters
// for bam file and reads
extern int min_flank_length;
extern int max_num_cigar;
extern int max_edit_distance;
extern int32_t min_bundle_gap;
extern int min_num_hits_in_bundle;
extern int min_num_splices_in_bundle;
extern uint32_t min_mapping_quality;
extern int32_t min_splice_boundary_hits;
extern bool uniquely_mapped_only;
extern bool use_second_alignment;

// for preview
extern bool preview_only;
extern int max_preview_reads;
extern int max_preview_spliced_reads;
extern int min_preview_spliced_reads;
extern int max_preview_umi_reads;
extern double preview_infer_ratio;
extern double insertsize_ave;
extern double insertsize_std;
extern int insertsize_median;
extern int insertsize_low;
extern int insertsize_high;
extern double insertsize_low_percentile;
extern double insertsize_high_percentile;
extern vector<double> insertsize_profile;

// for bridging
extern double min_bridging_score;
extern int max_num_path_nodes;
extern int dp_solution_size;
extern int dp_stack_size;
extern bool use_overlap_scoring;
extern int32_t max_clustering_flank; 
extern int32_t flank_tiny_length;
extern double flank_tiny_ratio;

// for identifying subgraphs
extern int32_t min_subregion_gap;
extern int32_t min_subregion_len;
extern int32_t min_subregion_max;
extern double min_subregion_ave;

// for amaranth
enum class amaranth_mode {REF, STAT_ONLY, MINI, ASSEMBLER};
enum class amaranth_strategy {LONG, SHORT, HEAVY};
enum class seq {SC, BULK, UNKNOWN};
extern amaranth_mode amaranthMode;
extern amaranth_strategy amaranthStrategy;
extern seq tech;
extern bool doesFastDnC;

// for subsetsum and router
extern int max_dp_table_size;
extern int min_router_count;

// for splice graph
extern double max_intron_contamination_coverage;
extern double min_surviving_edge_weight;
extern double max_decompose_error_ratio[7];
extern double min_transcript_numreads;
extern double min_transcript_coverage;
extern double min_guaranteed_edge_weight;

// for filtering transcripts
extern double min_single_exon_coverage;
extern double min_transcript_coverage_ratio; 
extern int min_transcript_length_base;
extern int min_transcript_length_increase;
extern int min_exon_length;
extern int max_num_exons;

// for simulation
extern int simulation_num_vertices;
extern int simulation_num_edges;
extern int simulation_max_edge_weight;

// input and output
extern string algo;
extern string input_file;
extern string ref_file;
extern string output_file;
extern string output_file1;

// umi & hybrid parameters
extern int min_umi_reads_bundle;
extern double min_umi_ratio_bundle;
extern bool both_umi_support;
extern int  min_umi_reads_start_exon;

// filtering & retention
extern bool remove_retained_intron;
extern int remove_dup;
extern bool use_filter;

// for controling
extern bool output_tex_files;
extern bool output_graphviz_files;
extern string fixed_gene_name;
extern int max_num_bundles;
extern int library_type;
extern int min_gtf_transcripts_num;
extern int batch_bundle_size;
extern int verbose;
extern int assemble_duplicates;
extern string version;

// parse arguments
int print_command_line(int argc, const char ** argv);
int parse_arguments(int argc, const char ** argv);
int print_parameters();
int print_copyright();
int print_help();

#endif