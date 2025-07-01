#include <bits/stdc++.h>
using namespace std;

/*============================================================================================================
Divide and Conquer Optimization

Description:
  • Optimizes DP of the form:
    DP(i, j) = min { DP(i-1, k-1) + Cost(k, j) } for k ≤ j
  • Requires the cost function to satisfy the "quadrangle inequality":
    F(a, c) + F(b, d) ≤ F(a, d) + F(b, c)  for a ≤ b ≤ c ≤ d
  • Ensures the monotonicity of the optimal k (OPT(i, j) ≤ OPT(i, j+1))
  • Reduces naive O(n·m²) DP complexity to approximately O(n·m·log(m)) using divide and conquer

How it works:
  - For each dp row, uses divide and conquer to efficiently find the optimal k for each position j
  - The search range for k narrows due to monotonicity, reducing computations

Parameters:
  - l, r: current range in dp_cur being computed
  - optl, optr: range of k candidates for dp transitions
  - dp_before: dp row i-1 (previous row)
  - dp_cur: dp row i (current row)

Usage:
  - Initialize dp_before with base cases
  - Call compute for each row to fill dp_cur efficiently
  - Assign dp_cur to dp_before for next iteration
============================================================================================================*/

int cost(int i, int j);

void compute(int l, int r, int optl, int optr, vector<int> &dp_before, vector<int> &dp_cur) {
    if (l > r) return;
    int mid = (l + r) >> 1;
    pair<int, int> best = {INT_MAX, -1};
    for (int k = optl; k <= min(mid, optr); k++) 
        best = min(best, {(k ? dp_before[k - 1] : 0) + cost(k, mid), k});
    dp_cur[mid] = best.first;
    int opt = best.second;
    compute(l, mid - 1, optl, opt, dp_before, dp_cur);    
    compute(mid + 1, r, opt, optr, dp_before, dp_cur);
}

void solve() {
    int n, m;
    vector<int> dp_before(m, 0), dp_cur(m, 0);
    for (int i = 0; i < m; i++) dp_before[i] = cost(0, i);
    for (int i = 1; i < n; i++) {
        compute(0,  m - 1, 0, m - 1, dp_before, dp_cur);
        dp_before = dp_cur;
    }
}

/*============================================================================================================
Knuth's Optimization

Description:
  • Optimizes DP of the form:
    DP(i, j) = min { DP(i, k) + DP(k + 1, j) + Cost(i, j) } for i ≤ k < j
  • Requires the cost function to satisfy the "quadrangle inequality" and the monotonicity condition:
    - Quadrangle Inequality: cost(a, c) + cost(b, d) ≤ cost(a, d) + cost(b, c)
    - Monotonicity: OPT(i, j - 1) ≤ OPT(i, j) ≤ OPT(i + 1, j)
  • Reduces time complexity from O(n³) to O(n²) by narrowing the search space for optimal k using the monotonicity property

Parameters:
  - n: Number of elements
  - cost(i, j): Function to compute the cost between indices i and j
  - dp: DP table to store the minimum cost for each subproblem
  - opt: Table to store the optimal k for each subproblem

Usage:
  - Initialize dp[i][i] and opt[i][i] for base cases
  - Iterate over increasing lengths of subarrays
  - For each subarray, compute dp[i][j] using the optimal k from opt[i][j-1] and opt[i+1][j]
  - The final answer will be stored in dp[0][n-1]
============================================================================================================*/

void solve_dp() {
    int n;
    vector<vector<int>> dp(n, vector<int>(n)), opt(n, vector<int>(n));
    for (int i = 0; i < n; i++) {
        opt[i][i] = i;
        // Initialize dp[i][i] according to the problem
    }
    for (int i = n - 2; i >= 0; i--) 
        for (int j = i + 1; j < n; j++) {
            int temp = INT_MAX, c = cost(i, j);
            for (int k = opt[i][j - 1]; k <= min(j - 1, opt[i + 1][j]); k++) 
                if (temp >= dp[i][k] + dp[k + 1][j] + c) {
                    opt[i][j] = k;
                    temp = dp[i][k] + dp[k + 1][j] + c;
                }
            dp[i][j] = temp;
        }
}