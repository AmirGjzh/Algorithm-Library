#include <bits/stdc++.h>
using namespace std;

/*--------------------------------------------------------------------------------------------------------
Divide and Conquer Optimization :
We can use this tecnique for DPs like this :
DP(i, j) = min {DP(i - 1, k - 1) + Cost(k, j)} for k <= j

We define OPT(i, j) = the k that minimize that expression
Also we assum that the cost function satisfies the "quadrangle inequality"
We can show that OPT(i, j) <= OPT(i, j + 1) also known as "monotonicity condition"
Now we use DC to solve our DP

dp_cur is the row that we are computing
dp_before is the last row that we computed
optl and optr, are decreasing the range we must go to find the best k
Our total Order goes O(n.m.log(m)) from O(n.m ^ 2)

A function satisfies the "quadrangle inequality" if and only if :
F(a, c) + F(b, d) <<<= F(a, d) + F(b, c)
By  A <<<= B, I mean A is more or equaly optimal than B
--------------------------------------------------------------------------------------------------------*/

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

/*--------------------------------------------------------------------------------------------------------
Knuth's Optimization :
We can use this tecnique for range DPs like this :
DP(i, j) = min {DP(i, k) + DP(k + 1, j)} + Cost(i, j)   for i <= k < j

We define OPT(i, j) = the k that minimize that expression
Also we assum that the cost function satisfies the "quadrangle inequality" and another condition
We can show that OPT(i, j - 1) <= OPT(i, j) <= OPT(i + 1, j)

Before we solve DP(i, j), we solve DP(i, j - 1) and DP(i + 1, j)
So then we know OPT(i, j - 1) <= OPT(i, j) <= OPT(i + 1, j)
We just for from OPT(i, j - 1) to OPT(i + 1, j)
Our total Order goes O(n ^ 3) from O(n ^ 2)

A function satisfies the "quadrangle inequality" if and only if :
a <= b <= c <= d
F(a, c) + F(b, d) <<<= F(a, d) + F(b, c)
By  A <<<= B, I mean A is more or equaly optimal than B

Another condition :
Cost(b, c) <= Cost(a, d)
--------------------------------------------------------------------------------------------------------*/

void solve_dp() {
    int n;
    vector<vector<int>> dp(n, vector<int>(n)), opt(n, vector<int>(n));
    for (int i = 0; i < n; i++) {
        opt[i][i] = i;
        // Initialize dp[i][i] according to the problem
    }
    for (int i = n - 2; i >= 0; i--) 
        for (int j = i + 1; j < n; j++) {
            int temp = INT_MAX;
            int c = cost(i, j);
            for (int k = opt[i][j - 1]; k <= min(j - 1, opt[i + 1][j]); k++) 
                if (temp >= dp[i][k] + dp[k + 1][j] + c) {
                    opt[i][j] = k;
                    temp = dp[i][k] + dp[k + 1][j] + c;
                }
            dp[i][j] = temp;
        }
}