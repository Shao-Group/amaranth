/*
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

#include "vertex_info.h"

vertex_info::vertex_info()
{
	stddev = 1.0;
	length = 0;
	sdist = -1;
	tdist = -1;
	type = -1;
	lpos = 0;
	rpos = 0;
	pos = 0;
	lstrand = '.';
	rstrand = '.';
	regional = false;
	umi_support = 0;
	cb_tags.clear();
}

vertex_info::vertex_info(const vertex_info &vi)
{
	stddev = vi.stddev;
	length = vi.length;
	sdist = vi.sdist;
	tdist = vi.tdist;
	type = vi.type;
	lpos = vi.lpos;
	rpos = vi.rpos;
	pos = vi.pos;
	lstrand = vi.lstrand;
	rstrand = vi.rstrand;
	regional = vi.regional;
	umi_support = vi.umi_support;
	cb_tags = vi.cb_tags;
}

int vertex_info::cell_support()
{
	return cb_tags.size();
}