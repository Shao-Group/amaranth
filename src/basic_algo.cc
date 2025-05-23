#include "basic_algo.h"

int basic_algo::smith_waterman(const vector<int>& s, const vector<int>& r, int gap_panelty, int mis_panelty, int match_reward)
{
    int m = s.size() + 1;
    int n = r.size() + 1;
    // initialization
	vector<vector<int> > opt;
	opt.resize(m, vector<int>(n));
    for (int i = 0; i < m; i++) 
    {
        for (int j = 0; j < n; j++) opt[i][j] = 0;
    }
    // body of dp table
    for (int i = 1; i < m; ++i)
    {
        for (int j = 1; j < n; ++j)
        {
			int cost = (s[i - 1] == r[j - 1]? match_reward : mis_panelty);
            int match    = opt[i - 1][j - 1] + cost;
            int gap_in_s = opt[i][j - 1]  + gap_panelty;
            int gap_in_r = opt[i - 1][j]  + gap_panelty;
            opt[i][j]  = maxv(vector<int> {0, match, gap_in_s, gap_in_r});
        }
    }
    return opt[m - 1][n - 1];
}

int basic_algo::needleman_wunsch(const vector<int>& s, const vector<int>& r,  int gap_panelty, int mis_panelty, int match_reward)
{
    int m = s.size() + 1;
    int n = r.size() + 1;
    // initialization
	vector<vector<int> > opt;
	opt.resize(m, vector<int>(n));
    for (int i = 0; i < m; i++) 
    {
        for (int j = 0; j < n; j++) opt[i][j] = 0;
    }
    for (int i = 0; i < m; i++) opt[i][0] = i * gap_panelty;
    for (int j = 0; j < n; j++) opt[0][j] = j * gap_panelty;
    // body of dp table
    for (int i = 1; i < m; ++i)
    {
        for (int j = 1; j < n; ++j)
        {
			int cost = (s[i - 1] == r[j - 1]? match_reward : mis_panelty);
            int match    = opt[i - 1][j - 1] + cost;
            int gap_in_s = opt[i][j - 1]  + gap_panelty;
            int gap_in_r = opt[i - 1][j]  + gap_panelty;
            opt[i][j]  = maxv(vector<int>{match, gap_in_s, gap_in_r});
        }
    }
    return opt[m - 1][n - 1];
}

/* TODO: test cases
*  align r to s, r should be all comsumed but s does not have to be comsumed up
*  e.g. r is read, s is genome; r is phasing path, s is complete path
*/ 

int basic_algo::ref_sw_query_nw(const vector<int>& ref_s, const vector<int>& query_r,  int gap_panelty, int mis_panelty, int match_reward)
{
    const vector<int>& s = ref_s;
    const vector<int>& r = query_r;
    int m = s.size() + 1;
    int n = r.size() + 1;
    // initialization
	vector<vector<int> > opt;
	opt.resize(m, vector<int>(n));
    for (int i = 0; i < m; i++) 
    {
        for (int j = 0; j < n; j++) opt[i][j] = 0;
    }
    for (int j = 0; j < n; j++) opt[0][j] = j * gap_panelty;
	
    // body of dp table
    for (int i = 1; i < m; ++i)
    {
        for (int j = 1; j < n; ++j)
        {
			int cost = (s[i - 1] == r[j - 1]? match_reward : mis_panelty);
            int match    = opt[i - 1][j - 1] + cost;
            int gap_in_s = opt[i][j - 1]  + gap_panelty;
            int gap_in_r = opt[i - 1][j]  + gap_panelty;
            opt[i][j]  = maxv(vector<int>{match, gap_in_s, gap_in_r});
        }
    }
    vector<int> optShortcutN(n);
    for(int i = 0; i < n; i ++) optShortcutN[i] = opt[m - 1][i];
    int optN = maxv(optShortcutN);
    return optN;
}