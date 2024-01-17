#include "basic_algo.h"

int basic_algo::smith_waterman(vector<int>& s, vector<int>& r)
{
    int m = s.size() + 1;
    int n = r.size() + 1;
	const int GAP_PANELTY = -1;
	const int MIS_PANELTY = -1;
    const int MATCH_REWARD = 1;

    // initialization
	vector<vector<int> > opt;
	opt.resize(m, vector<int>(n));
    for (int i = 0; i < m; i++) for (int j = 0; j < n; j++) opt[i][j] = 0;
    for (int i = 0; i < m; i++) opt[i][0] = i * GAP_PANELTY;
    for (int j = 0; j < n; j++) opt[0][j] = j * GAP_PANELTY;
	
    // body of dp table
    for (int i = 1; i < m; ++i)
    {
        for (int j = 1; j < n; ++j)
        {
			int cost = (s[i - 1] == r[j - 1]? MIS_PANELTY : MIS_PANELTY);
            int match    = opt[i - 1][j - 1] + cost;
            int gap_in_s = opt[i][j - 1]  + GAP_PANELTY;
            int gap_in_r = opt[i - 1][j]  + GAP_PANELTY;
            opt[i][j]  = maxv(vector<int> {0, match, gap_in_s, gap_in_r});
        }
    }
    
    return opt[m - 1][n - 1];
}

int basic_algo::needleman_wunsch(vector<int>& s, vector<int>& r)
{
    int m = s.size() + 1;
    int n = r.size() + 1;
	const int GAP_PANELTY = -1;
	const int MIS_PANELTY = -1;
    const int MATCH_REWARD = 1;

    // initialization
	vector<vector<int> > opt;
	opt.resize(m, vector<int>(n));
    for (int i = 0; i < m; i++) for (int j = 0; j < n; j++) opt[i][j] = 0;
    for (int i = 0; i < m; i++) opt[i][0] = i * GAP_PANELTY;
    for (int j = 0; j < n; j++) opt[0][j] = j * GAP_PANELTY;
	
    // body of dp table
    for (int i = 1; i < m; ++i)
    {
        for (int j = 1; j < n; ++j)
        {
			int cost = (s[i - 1] == r[j - 1]? MIS_PANELTY : MIS_PANELTY);
            int match    = opt[i - 1][j - 1] + cost;
            int gap_in_s = opt[i][j - 1]  + GAP_PANELTY;
            int gap_in_r = opt[i - 1][j]  + GAP_PANELTY;
            // opt[i][j]  = max(match, max(gap_in_s, gap_in_r));
            opt[i][j]  = maxv(vector<int>{match, gap_in_s, gap_in_r});
        }
    }
    
    return opt[m - 1][n - 1];
}
