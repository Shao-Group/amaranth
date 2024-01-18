
#ifndef __BASIC_ALGO_H__
#define __BASIC_ALGO_H__

#include<vector>
#include<cassert>

using namespace::std;

namespace basic_algo
{
    int needleman_wunsch(vector<int>& s, vector<int>& r, int gap = -1, int mis_panelty = -1, int match = 1);
    int smith_waterman(vector<int>& s, vector<int>& r,  int gap = -1, int mis_panelty = -1, int match = 1);    
}

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