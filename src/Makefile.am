bin_PROGRAMS = amaranth


GTF_INCLUDE = $(top_srcdir)/lib/gtf
UTIL_INCLUDE = $(top_srcdir)/lib/util
GRAPH_INCLUDE = $(top_srcdir)/lib/graph

GTF_LIB = $(top_builddir)/lib/gtf
UTIL_LIB = $(top_builddir)/lib/util
GRAPH_LIB = $(top_builddir)/lib/graph

amaranth_CPPFLAGS = -std=c++20 -I$(GTF_INCLUDE) -I$(GRAPH_INCLUDE) -I$(UTIL_INCLUDE)
amaranth_LDFLAGS = -pthread -L$(GTF_LIB) -L$(GRAPH_LIB) -L$(UTIL_LIB)
amaranth_LDADD = -lgtf -lgraph -lutil

amaranth_SOURCES = splice_graph.h splice_graph.cc \
				   super_graph.h super_graph.cc \
				   sgraph_compare.h sgraph_compare.cc \
				   vertex_info.h vertex_info.cc \
				   edge_info.h edge_info.cc \
				   interval_map.h interval_map.cc \
				   config.h config.cc \
				   hit.h hit.cc \
				   partial_exon.h partial_exon.cc \
				   hyper_set.h hyper_set.cc \
				   subsetsum.h subsetsum.cc \
				   router.h router.cc \
				   region.h region.cc \
				   junction.h junction.cc \
				   fragment.h fragment.cc \
				   bundle_base.h bundle_base.cc \
				   bundle_bridge.h bundle_bridge.cc \
				   bridger.h bridger.cc \
				   bundle.h bundle.cc \
				   fcluster.h fcluster.cc \
				   path.h path.cc \
				   equation.h equation.cc \
				   gtf.h gtf.cc \
				   previewer.h previewer.cc \
				   assembler.h assembler.cc \
				   transcript_set.h transcript_set.cc \
				   filter.h filter.cc \
				   amaranth_aux.hpp \
				   amaranth.h amaranth.cc \
				   scallop.h scallop.cc \
				   analyzer.h analyzer.cc \
				   basic_algo.h basic_algo.cc \
				   main.cc
