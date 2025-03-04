/*
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

#include <cassert>
#include <cstdio>
#include <cmath>
#include <climits>

#include "bundle_base.h"

bundle_base::bundle_base()
{
	tid = -1;
	chrm = "";
	lpos = 1 << 30;
	rpos = 0;
	strand = '.';
}

bundle_base::~bundle_base()
{}

int bundle_base::add_hit(const hit &ht)
{
	// store new hit
	hits.push_back(ht);

	// calcuate the boundaries on reference
	if(ht.pos < lpos) lpos = ht.pos;
	if(ht.rpos > rpos) rpos = ht.rpos;

	// try to include more paired-end reads
	int32_t p = ht.rpos;
	if(ht.mpos > ht.rpos && ht.mpos <= ht.rpos + 100000) p = ht.mpos;
	if(p > rpos) rpos = p;

	// set tid
	if(tid == -1) tid = ht.tid;
	assert(tid == ht.tid);

	// set strand
	if(hits.size() <= 1) strand = ht.strand;
	assert(strand == ht.strand);

	// DEBUG
	/*
	if(strand != ht.strand)
	{
		printf("strand = %c, ht.strand = %c, ht.xs = %c,\n", strand, ht.strand, ht.xs);
	}
	*/
	return 0;
}

// remove duplicated hits
// duplicated hits are defined as hit-pairs (wrt qname) with the same alignment and CIGAR
int bundle_base::rm_duplicated_reads()
{
	// vector <qname, hit-index>
	vector<pair<string, int>> qh;	
	for(int i = 0; i < hits.size(); i++)
	{	
		qh.push_back(make_pair(hits[i].qname, i));
	}
	sort(qh.begin(), qh.end());

	// map <alignment x CIGAR -> qname>
	map<size_t, string> m;
	string current_qname = "";
	stringstream current_info;
	current_info.str("");
	current_info.clear();
	for(const auto& [qname, i]: qh)
	{
		if(qname != current_qname && current_info.str() != "")
		{
			int hh = string_hash(current_info.str());
			if(m.find(hh) == m.end())
			{
				// printf("current_info = %s, hh = %d, added\n", current_info.str().c_str(), hh); // debug
				m.insert(make_pair(hh, current_qname));
			}
			// else
			// {
			// 	printf("current_info = %s, hh = %d, not added\n", current_info.str().c_str(), hh); // debug
			// }
			current_qname = qname;
			current_info.str("");
			current_info.clear();
		}
		// get all hits with the same qname and hash their alignment and CIGAR
		current_info << hits[i].pos  << "-" << hits[i].cigar_str << ".";
	}
	int hh = string_hash(current_info.str());
	if(m.find(hh) == m.end() && current_qname != "")  
	{
		m.insert(make_pair(hh, current_qname));
	}
	
	// sort qnames, keep survived ones
	vector<string> qkeep;
	for(const auto& [hh, qname]: m)	qkeep.push_back(qname);
	sort(qkeep.begin(), qkeep.end());

	// keep hits with survived qnames
	vector<int> ikeep;
	int p1 = 0, p2 = 0;
	while(p1 < qh.size() && p2 < qkeep.size())
	{
		if(qh[p1].first == qkeep[p2])
		{
			ikeep.push_back(qh[p1].second);
			p1++;
		}
		else if(qh[p1].first < qkeep[p2]) p1++;
		else p2++;
	}
	sort(ikeep.begin(), ikeep.end());
	vector<hit> v;
	for (int i = 0; i < ikeep.size(); i++) 
	{
		if (i >= 1) assert(ikeep[i] > ikeep[i -1]);
		assert(ikeep[i] < hits.size());
		v.push_back(hits[ikeep[i]]);
	}
	hits = v;
	return 0;
}

int bundle_base::build_maps()
{
	for (const auto& ht: hits)
	{

		for(int k = 0; k < ht.itvm.size(); k++)
		{
			int32_t s = high32(ht.itvm[k]);
			int32_t t = low32(ht.itvm[k]);
			//printf(" add interval %d-%d\n", s, t);
			mmap += make_pair(ROI(s, t), 1);
		}

		for(int k = 0; k < ht.itvi.size(); k++)
		{
			int32_t s = high32(ht.itvi[k]);
			int32_t t = low32(ht.itvi[k]);
			imap += make_pair(ROI(s, t), 1);
		}

		for(int k = 0; k < ht.itvd.size(); k++)
		{
			int32_t s = high32(ht.itvd[k]);
			int32_t t = low32(ht.itvd[k]);
			imap += make_pair(ROI(s, t), 1);
		}
	}
	return 0;
}

bool bundle_base::overlap(const hit &ht) const
{
	if(mmap.find(ROI(ht.pos, ht.pos + 1)) != mmap.end()) return true;
	if(mmap.find(ROI(ht.rpos - 1, ht.rpos)) != mmap.end()) return true;
	return false;
}

int bundle_base::clear()
{
	tid = -1;
	chrm = "";
	lpos = 1 << 30;
	rpos = 0;
	strand = '.';
	hits.clear();
	mmap.clear();
	imap.clear();
	return 0;
}

