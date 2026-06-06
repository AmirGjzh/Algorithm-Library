#include <bits/stdc++.h>
using namespace std;
using ll = long long int;

/*============================================================================================================
Classic Dynamic Programming & Sequence Algorithms

Contents:
  • Knapsack family (0–1, Unbounded, Multiple)
  • Longest Increasing Subsequence (LIS)
  • Longest Common Subsequence (LCS)
  • Largest Zero Submatrix in a binary grid

Knapsack:
  • Goal: maximize total value under weight capacity W
  • 0–1 Knapsack:
    – Each item used once
    – DP over capacity, iterating weights backward
    – Time O(n·W), Space O(W)
  • Unbounded Knapsack:
    – Unlimited use per item
    – DP over capacity, iterating forward
    – Time O(n·W)
  • Multiple Knapsack:
    – Bounded count per item
    – Binary decomposition converts it to 0–1 knapsack
    – Effective item count ≈ Σ log(k_i)

LIS (Longest Increasing Subsequence):
  • Finds longest strictly increasing subsequence (not contiguous)
  • O(n²) DP version:
    – dp[i]: LIS ending at i
    – Uses parent array to reconstruct sequence
  • O(n·log(n)) optimized version:
    – dp[len]: smallest possible tail of an LIS of length len
    – Uses binary search and parent tracking for reconstruction
  • Variants:
    – Non-decreasing LIS
    – Counting number of LIS
    – Minimum number of non-increasing subsequences covering sequence = LIS length

LCS (Longest Common Subsequence):
  • Finds longest subsequence common to two strings
  • dp[i][j]: LCS length of prefixes s[0..i-1], t[0..j-1]
  • Classic 2D DP with reconstruction
  • Time O(n·m), Space O(n·m)

Largest Zero Submatrix:
  • Finds the largest rectangle consisting only of zeros in a binary matrix
  • Treats each row as the bottom of a histogram
  • d[j]: last row index containing a 1 in column j
  • Uses monotonic stacks to find max rectangle per row
  • Time O(n·m)

Notes:
  • All implementations are standard, static, and contest-safe
  • Focus is on clarity and correctness rather than micro-optimizations
============================================================================================================*/

ll knapsack_01(const int W, const vector<int> &w, const vector<int> &v) {
    const int n = static_cast<int>(w.size());
    vector<ll> dp(W + 1, 0);
    for (int i = 0; i < n; i++)
        for (int j = W; j >= w[i]; j--)
            dp[j] = max(dp[j], dp[j - w[i]] + v[i]);
    return dp[W];
}

ll unbounded_knapsack(const int W, const vector<int> &w, const vector<int> &v) {
    const int n = static_cast<int>(w.size());
    vector<ll> dp(W + 1, 0);
    for (int i = 0; i < n; i++)
        for (int j = w[i]; j <= W; j++)
            dp[j] = max(dp[j], dp[j - w[i]] + v[i]);
    return dp[W];
}

ll multiple_knapsack(const int W, const vector<int> &w, const vector<int> &v, vector<int> k) {
    const int n = static_cast<int>(w.size());
    vector<int> new_w, new_v;
    for (int i = 0; i < n; i++)
        for (int p = 1; k[i] > 0; p <<= 1) {
            const int take = min(p, k[i]);
            new_w.push_back(take * w[i]);
            new_v.push_back(take * v[i]);
            k[i] -= take;
        }
    return knapsack_01(W, new_w, new_v);
}

vector<ll> LIS_bad(const vector<ll> &a) {
    const int n = static_cast<int>(a.size());
    int ans = 1, pos = 0;
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
    const int n = static_cast<int>(a.size());
    constexpr ll INF = LLONG_MAX;
    int pos = -1;
    vector dp(n + 1, INF);
    vector<int> ind(n + 1, -1), prev(n, -1);
    dp[0] = -INF;
    for (int i = 0; i < n; i++) {
        if (const int l = upper_bound(dp.begin(), dp.end(), a[i]) - dp.begin(); dp[l - 1] < a[i] and a[i] < dp[l]) {
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

string LCS(const string &s, const string &t) {
    string lcs;
    const int n = static_cast<int>(s.size()), m = static_cast<int>(t.size());
    vector dp(n + 1, vector(m + 1, 0));
    for (int i = 1; i <= n; i++)
        for (int j = 1; j <= m; j++)
            if (s[i - 1] == t[j - 1]) dp[i][j] = dp[i - 1][j - 1] + 1;
            else dp[i][j] = max(dp[i - 1][j], dp[i][j - 1]);
    int i = n, j = m;
    while (i > 0 and j > 0)
        if (s[i - 1] == t[j - 1]) {
            lcs.insert(lcs.begin(), s[i - 1]);
            i--, j--;
        }
        else if (dp[i - 1][j] > dp[i][j - 1]) i--;
        else j--;
    return lcs;
}

int zero_submatrix(const vector<vector<int>> &a) {
    const int n = static_cast<int>(a.size()), m = static_cast<int>(a[0].size());
    int ans = 0;
    vector<int> d(m, -1), l(m), r(m);
    stack<int> st;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++)
            if (a[i][j] == 1) d[j] = i;
        while (!st.empty()) st.pop();
        for (int j = 0; j < m; j++) {
            while (!st.empty() and d[st.top()] <= d[j]) st.pop();
            l[j] = st.empty() ? -1 : st.top();
            st.push(j);
        }
        while (!st.empty()) st.pop();
        for (int j = m - 1; j >= 0; j--) {
            while (!st.empty() and d[st.top()] <= d[j]) st.pop();
            r[j] = st.empty() ? m : st.top();
            st.push(j);
        }
        for (int j = 0; j < m; j++) ans = max(ans, (i - d[j]) * (r[j] - l[j] - 1));
    }
    return ans;
}
