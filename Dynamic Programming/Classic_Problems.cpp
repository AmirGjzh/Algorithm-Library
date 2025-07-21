#include <bits/stdc++.h>
using namespace std;
using ll = long long int;

/*============================================================================================================
Knapsack Problem

Description:
  • Given items with weight wi and value vi, and a capacity W, choose a subset to maximize total value 
    without exceeding W

Variants:
  1. 0–1 Knapsack (each item used at most once)
    – dp(i, j) = max value from first i items with capacity j:
      dp[i][j] = max(dp[i-1][j], dp[i-1][j - wi] + vi)
    – Time: O(n·W), Space: O(n·W)
  2. Unbounded Knapsack (unlimited uses per item)
    – Transition: dp[i][j] = max(dp[i-1][j], dp[i][j - wi] + vi)
    – Time: O(n·W)
  3. Multiple Knapsack (bounded count ki per item)
    – Use binary-splitting: decompose ki into sums of powers of two → convert to 0–1 knapsack
    – Effective item-count n′ ≈ Σ log ki; Time: O(n′·W)

Applications:
  • Resource allocation, scheduling, budgeting, portfolio optimization
============================================================================================================*/

ll knapsack_01(int W, const vector<int> &w, const vector<int> &v) {
    int n = w.size();
    vector<ll> dp(W + 1, 0);
    for (int i = 0; i < n; i++) 
        for (int j = W; j >= w[i]; j--)
            dp[j] = max(dp[j], dp[j - w[i]] + v[i]);
    return dp[W];
}

ll unbounded_knapsack(int W, const vector<int> &w, const vector<int> &v) {
    int n = w.size();
    vector<ll> dp(W + 1, 0);
    for (int i = 0; i < n; i++) 
        for (int j = w[i]; j <= W; j++) 
            dp[j] = max(dp[j], dp[j - w[i]] + v[i]);
    return dp[W];    
}

ll multiple_knapsack(int W, const vector<int> &w, const vector<int> &v, vector<int> k) {
    int n = w.size();
    vector<int> new_w, new_v;
    for (int i = 0; i < n; i++) 
        for (int p = 1; k[i] > 0; p <<= 1) {
            int take = min(p, k[i]);
            new_w.push_back(take * w[i]);
            new_v.push_back(take * v[i]);
            k[i] -= take;
        }
    return knapsack_01(W, new_w, new_v);
}

/*============================================================================================================
Longest Increasing Subsequence (LIS)

Description:
  • Given a sequence, find the longest strictly increasing subsequence (not necessarily contiguous)

Methods:
  1. DP O(n²):
    – dp[i]: length of LIS ending at index i
    – Transition: dp[i] = max(dp[j]) + 1 for all j < i with a[j] < a[i]
  2. DP + Binary Search O(n.log(n)):
    – dp[len]: smallest tail value of an increasing subsequence of length len
    – For each element, binary search dp to update subsequence length efficiently
  3. Data Structures (Segment Tree / Fenwick Tree):
    – Use DS to optimize max queries in DP, useful for advanced variants

Variants & Tasks:
  • Longest non-decreasing subsequence: use ≤ instead of <, adjust binary search
  • Counting number of LIS: DP with counting, handle equal-length subsequences
  • Minimum number of non-increasing subsequences covering sequence: equals LIS length
    – The minimal number of colors needed so that each colored group is non-increasing
    – Answer = size of LIS
    – Greedy coloring is possible

Applications:
  • Pattern recognition, data analysis, bioinformatics, stock market trend analysis
============================================================================================================*/

vector<ll> LIS_bad(const vector<ll> &a) {
    int n = a.size(), ans = 1, pos = 0;
    vector<int> dp(n, 1), prev(n, -1);
    for (int i = 0; i < n; i++) 
        for (int j = 0; j < i; j++) 
            if (a[j] < a[i] and dp[i] < dp[j] + 1) 
                dp[i] = dp[j] + 1, prev[i] = j;
    for (int i = 1; i < n; i++) 
        if (dp[i] > ans) 
            ans = dp[i], pos = i;
    vector<ll> lis;
    while (pos != -1) 
        lis.push_back(a[pos]), pos = prev[pos];
    reverse(lis.begin(), lis.end());
    return lis;
}

vector<ll> LIS(const vector<ll> &a) {
    int n = a.size(), pos = -1, INF = LLONG_MAX;
    vector<ll> dp(n + 1, INF);
    vector<int> ind(n + 1, -1), prev(n, -1);
    dp[0] = -INF;
    for (int i = 0; i < n; i++) {
        int l = upper_bound(dp.begin(), dp.end(), a[i]) - dp.begin();
        if (dp[l - 1] < a[i] and a[i] < dp[l]) {
            dp[l] = a[i];
            ind[l] = i;
            prev[i] = ind[l - 1];
        }
    }
    for (int i = 0; i <= n; i++) 
        if (dp[i] < INF) pos = ind[i];
    vector<ll> lis;
    while (pos != -1) lis.push_back(a[pos]), pos = prev[pos];
    reverse(lis.begin(), lis.end());
    return lis;
}

/*============================================================================================================
Longest Common Subsequence (LCS)

Description:
  • Given two strings s and t, find the longest subsequence common to both
  • A subsequence is a sequence that appears in the same order, not necessarily contiguous

DP Definition:
  • dp[i][j] = length of LCS of s[0..i-1] and t[0..j-1]

Transition:
  • If s[i-1] == t[j-1], dp[i][j] = dp[i-1][j-1] + 1
  • Else dp[i][j] = max(dp[i-1][j], dp[i][j-1])

Time Complexity:
  • O(n.m), where n = length of s, m = length of t

Extensions:
  • LCS of k strings can be found with k-dimensional DP, but complexity grows exponentially

Applications:
  • Bioinformatics (DNA sequence alignment), text comparison, diff tools, version control systems
============================================================================================================*/

string LCS(const string &s, const string &t) {
    string lcs = "";
    int n = s.size(), m = t.size();
    vector<vector<int>> dp(n + 1, vector<int>(m + 1, 0));
    for (int i = 1; i <= n; i++) 
        for (int j = 1; j <= m; j++) 
            if (s[i - 1] == t[j - 1]) dp[i][j] = dp[i - 1][j - 1] + 1;
            else dp[i][j] = max(dp[i - 1][j], dp[i][j - 1]);
    int i = n, j = m;
    while (i > 0 and j > 0) 
        if (s[i - 1] == t[j - 1]) {
            lcs = s[i - 1] + lcs;
            i--, j--;
        }
        else if (dp[i - 1][j] > dp[i][j - 1]) i--;
        else j--;
    return lcs; 
}

/*============================================================================================================
Largest Zero Submatrix

Description:
  • Given a binary matrix, find the largest rectangular submatrix consisting entirely of zeros
  • For each row fixed as the bottom boundary, determine the heights of continuous zero-columns above
  • Use a stack-based approach to compute the largest rectangle area for each row's histogram representation

Key variables:
  • d[j]: Nearest row above (or equal) where column j has a '1' (or -1 if none)
  • l[j], r[j]: Indices for left and right boundaries of the maximal rectangle for column j

Algorithm steps:
  1. For each row i as the bottom boundary:
    - Update d[] for columns where '1' is present
  2. Use a stack to find previous smaller elements for l[] and next smaller for r[]
  3. Calculate area for each column as height * width = (i - d[j]) * (r[j] - l[j] - 1)
  4. Keep track of the maximum area

Time Complexity:
  • O(n.m), where n and m are dimensions of the matrix

Applications:
  • Image processing, free space detection in grids, pattern recognition
============================================================================================================*/

int zero_submatrix(const vector<vector<int>> &a) {
    int n = a.size(), m = a[0].size(), ans = 0;
    vector<int> d(m, -1), l(m), r(m);
    stack<int> st;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) 
            if (a[i][j] == 1) d[j] = i;
        while (st.size()) st.pop();    
        for (int j = 0; j < m; j++) {
            while (st.size() and d[st.top()] <= d[j]) st.pop();
            l[j] = st.empty() ? -1 : st.top();
            st.push(j);    
        }
        while (st.size()) st.pop();
        for (int j = m - 1; j >= 0; j--) {
            while (st.size() and d[st.top()] <= d[j]) st.pop();
            r[j] = st.empty() ? m : st.top();
            st.push(j); 
        }
        for (int j = 0; j < m; j++) ans = max(ans, (i - d[j]) * (r[j] - l[j] - 1));
    }
    return ans;
}
