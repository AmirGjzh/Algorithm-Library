#include <bits/stdc++.h>
using namespace std;
using ll = long long int;

/*============================================================================================================
DP Optimizations: Divide & Conquer Optimization + Knuth Optimization

Purpose:
  • Speed up specific DP transitions by shrinking the search space for the optimal split point k
  • Both rely on monotonicity of optimal decisions, which comes from properties of the cost function

1) Divide and Conquer Optimization

Applies to DP of the form:
  • DP(i, j) = min over k ≤ j { DP(i-1, k-1) + cost(k, j) }

Requirements:
  • Cost function satisfies quadrangle inequality
  • Optimal split positions are monotonic: OPT(i, j) ≤ OPT(i, j+1)

Idea:
  • Compute each DP row using divide & conquer
  • For a segment [l, r], search k only in [optl, optr]
  • Best k at mid limits the search range for left/right halves

Complexity:
  • From O(n·m²) → O(n·m·log(m))

Notes:
  • dp_before = DP row i-1
  • dp_cur     = DP row i
  • cost(k, j) must be fast (usually O(1))

2) Knuth Optimization

Applies to interval DP of the form:
  • DP(i, j) = min over i ≤ k < j { DP(i, k) + DP(k+1, j) + cost(i, j) }

Requirements:
  • Quadrangle inequality on cost
  • Strong monotonicity: OPT(i, j-1) ≤ OPT(i, j) ≤ OPT(i+1, j)

Idea:
  • When computing dp[i][j], only try k in: [ opt[i][j-1], opt[i+1][j] ]
  • Store optimal k to reuse later

Complexity:
  • From O(n³) → O(n²)

Notes:
  • dp[i][i] are base cases
  • opt[i][j] stores the best split point
  • Common in optimal BST, merging intervals, parsing problems
============================================================================================================*/

ll cost(int i, int j);

void compute(const int l, const int r, const int optl, const int optr, const vector<ll> &dp_before, vector<ll> &dp_cur) {
    if (l > r) return;
    const int mid = (l + r) >> 1;
    pair best = {LLONG_MAX, -1};
    for (int k = optl; k <= min(mid, optr); k++) 
        best = min(best, {(k ? dp_before[k - 1] : 0) + cost(k, mid), k});
    dp_cur[mid] = best.first;
    const int opt = best.second;
    compute(l, mid - 1, optl, opt, dp_before, dp_cur);    
    compute(mid + 1, r, opt, optr, dp_before, dp_cur);
}

void solve() {
    int n, m;
    vector<ll> dp_before(m), dp_cur(m);
    for (int i = 0; i < m; i++) dp_before[i] = cost(0, i);
    for (int i = 1; i < n; i++) {
        compute(0,  m - 1, 0, m - 1, dp_before, dp_cur);
        dp_before.swap(dp_cur);
    }
}

void solve_dp() {
    int n;
    vector dp(n, vector<ll>(n));
    vector opt(n, vector<int>(n));
    for (int i = 0; i < n; i++) {
        opt[i][i] = i;
        // Initialize dp[i][i] according to the problem
    }
    for (int i = n - 2; i >= 0; i--) 
        for (int j = i + 1; j < n; j++) {
            ll temp = LLONG_MAX;
            const ll c = cost(i, j);
            for (int k = opt[i][j - 1]; k <= min(j - 1, opt[i + 1][j]); k++) 
                if (temp >= dp[i][k] + dp[k + 1][j] + c) {
                    opt[i][j] = k;
                    temp = dp[i][k] + dp[k + 1][j] + c;
                }
            dp[i][j] = temp;
        }
}
