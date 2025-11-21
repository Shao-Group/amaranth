/*
Part of Amaranth Transcript Assembler
(c) 2025 by Xiaofei Carl Zang, Mingfu Shao, and The Pennsylvania State University.
See LICENSE for licensing.
*/

#ifndef __BASIC_ALGO_H__
#define __BASIC_ALGO_H__

#include<vector>
#include<cassert>

using namespace::std;

namespace basic_algo
{
    int needleman_wunsch(const vector<int>& s, const vector<int>& r, int gap = -1, int mis = -1, int match = 1);
    int smith_waterman(const vector<int>& complate_s, const vector<int>& partial_r,  int gap = -1, int mis = -1, int match = 1);    
    int ref_sw_query_nw(const vector<int>& s, const vector<int>& r, int gap = -1, int mis = -1, int match = 1);
}

// max value of a vector
template<typename T>
T maxv(vector<T> vt)
{
    assert(vt.size() > 0);
    T maxt = vt.front();
    for(const auto& t: vt)
    {
        if(maxt >= t) continue;
        maxt = t;
    }
    return maxt;
}

#endif