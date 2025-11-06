/*
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

#ifndef __VERTEX_INFO__
#define __VERTEX_INFO__

#include <stdint.h>
#include <set>
#include <string>
#include <map>

using namespace std;

class vertex_info
{
public:
	vertex_info();
	vertex_info(const vertex_info &vi);

public:
	int32_t pos;		// position
	int32_t lpos;		// left position
	int32_t rpos;		// right position
	double stddev;		// standard deviation of read coverage
	int length;			// length of this partial exon
	int sdist;			// shortest distance to s
	int tdist;			// shortest distance to t
	int type;			// for various usage
	char lstrand;		// left side strand
	char rstrand;		// right side strand	
	bool regional;		// if a vertex is regional
	int umi_support;
	map<string, set<string>> cb_tags;  // <CB, <UMI> >

	int cell_support(); // return number cell support (size of cb_tags)
};

#endif
