#include <bits/stdc++.h>
using namespace std;

/*--------------------------------------------------------------------------------------------------------
Knapsack problem :
We have a knapsack of maximum weight W, and n items that the i'th item has value vi and weight wi
We want to take some of the itmes in a way that we reach the maximum value

0-1 :
Difinition : dp(i, j) = maximum value we can reach from the i first items, and maximum wieght j
Update : dp(i, j) = max(dp(i - 1, j), dp(i - 1, j - wi) + vi)
Order = O(n.W)

Unbounded :
Definition : dp(i, j) = maximum value we can reach from the i first items, and maximum wieght j
Update : dp(i, j) = max(dp(i - 1, j), dp(i, j - wi) + vi)
Order = O(n.W)

Multiple :
This is just a 0-1 version, we must make it look like that. We know that the i'th element has ki of itself,
we can use binary grouping to convert "ki items of i'th type" to "Log(ki) items of groups of i'th type", 
in this way we can collect any number of i'th type, by collecting some of that Log(ki) groups
Difinition : dp(i, j) = maximum value we can reach from the i first items, and maximum wieght j
Update : dp(i, j) = max(dp(i - 1, j), dp(i - 1, j - wi) + vi)
Order = O(n'.W) here n' is different from n (Log(k1) + Log(k2) + ... + Log(k2))
--------------------------------------------------------------------------------------------------------*/

int knapsack_01(int W, vector<int> &w, vector<int> &v) {
    int n = w.size();
    vector<vector<int>> dp(n + 1, vector<int>(W + 1, 0));
    for (int i = 1; i <= n; i++) 
        for (int j = 0; j <= W; j++) {
            dp[i][j] = dp[i - 1][j];
            if (j >= w[i - 1])
                dp[i][j] = max(dp[i][j], dp[i - 1][j - w[i - 1]] + v[i - 1]);
        }
    return dp[n][W];
}

int unbounded_knapsack(int W, vector<int> &w, vector<int> &v) {
    int n = w.size();
    vector<vector<int>> dp(n + 1, vector<int>(W + 1, 0));
    for (int i = 1; i <= n; i++) 
        for (int j = 0; j <= W; j++) {
            dp[i][j] = dp[i - 1][j];
            if (j >= w[i - 1])
                dp[i][j] = max(dp[i][j], dp[i][j - w[i - 1]] + v[i - 1]);
        }
    return dp[n][W];    
}

int multiple_knapsack(int W, vector<int> &w, vector<int> &v, vector<int> &k) {
    int n = w.size();
    vector<int> new_w, new_v;
    for (int i = 0; i < n; i++) {
        int c = 1;
        while (k[i] > c) {
            k[i] -= c;
            new_w.push_back(c * w[i]);
            new_v.push_back(c * v[i]);
            c <<= 1;
        }
        new_w.push_back(k[i] * w[i]);
        new_v.push_back(k[i] * v[i]);
    }
    return knapsack_01(W, new_w, new_v);
}

/*--------------------------------------------------------------------------------------------------------
Longest Increasing Subsequence :
First way :
DP in O(n ^ 2)
Definition : dp(i) = size of LIS that ends in i'th place
Update : "dp(i) = max(dp(j)) + 1" that j < i and a[j] < a[i]

Second way :
DP and BS with O(n.Log(n))
Definition : dp(i) = smallest element that, a subsequence of size i ends with (INF if there is no answer)
Update : we iterate over our array and update the DP

Third way :
using DSs like segtree or fenwick
we go with LIS_bad, but for finding the max a[j], we use that DS

Some tasks :
1 - Longest non-decreasing subsequence :
It's the same as before, but we have to use <= instead of <
and in binary search part we have to be carefull

2 - Number of longest increasing subsequences
First way => along side of the answer to dp(i), we save the count of the answer too
it meanse we have to check the == situations
Second way => not possible
Third way => finding the max  and its frequency with seg tree :)

3 - Smallest number of non-increasing subsequences covering a sequence
Definition : smallset number of colors need to color an array, so that every group of the same color are non-increasing
The answer is the size of LIS :)
for coloring, we can greedily color the array
--------------------------------------------------------------------------------------------------------*/

vector<int> LIS_bad(vector<int> &a) {
    int n = a.size(), ans = 1, pos = 0;
    vector<int> dp(n, 1), prev(n, -1);
    for (int i = 0; i < n; i++) 
        for (int j = 0; j < i; j++) 
            if (a[j] < a[i] and dp[i] < dp[j] + 1) 
                dp[i] = dp[j] + 1, prev[i] = j;
    for (int i = 1; i < n; i++) 
        if (dp[i] > ans) 
            ans = dp[i], pos = i;
    vector<int> lis;
    while (pos != -1) 
        lis.push_back(a[pos]), pos = prev[pos];
    reverse(lis.begin(), lis.end());
    return lis;
}

vector<int> LIS(vector<int> &a) {
    int n = a.size(), pos = -1, INF = 1e9 + 10;
    vector<int> dp(n + 1, INF), ind(n + 1, -1), prev(n, -1);
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
        if (dp[i] < INF) 
            pos = ind[i];
    vector<int> lis;
    while (pos != -1) 
        lis.push_back(a[pos]), pos = prev[pos];
    reverse(lis.begin(), lis.end());
    return lis;
}

/*--------------------------------------------------------------------------------------------------------
Longest Common Subsequence :
Given 2 strings, find the LCS of them in O(n ^ 2)
Definition : dp(i, j) = LCS of s[0:i] and t[0:j]
Update : dp(i, j) = if s[i] == t[j] then dp(i - 1, j - 1) + 1 else max(dp(i - 1, j), dp(i, j - 1))
You can find LCS of k strings with a k-dimantion DP
--------------------------------------------------------------------------------------------------------*/

string LCS(string &s, string &t) {
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

/*--------------------------------------------------------------------------------------------------------
Finding the largest zero submatrix :
We fix the bottom of the submatrix and use the stack (like the classic problem) to find the answer
d[i][j] => the nearest row above a[i][j] that is 1 (it can be i)
--------------------------------------------------------------------------------------------------------*/

int zero_submatrix(vector<vector<int>> &a) {
    int n = a.size(), m = a[0].size(), ans = 0;
    vector<int> d(m, -1), l(m), r(m);
    stack<int> st;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) 
            if (a[i][j] == 1) d[j] = i;
        for (int j = 0; j < m; j++) {
            while (st.size() and d[st.top()] <= d[j]) 
                st.pop();
            l[j] = st.empty() ? -1 : st.top();
            st.push(j);    
        }
        while (st.size()) st.pop();
        for (int j = m - 1; j >= 0; j--) {
            while (st.size() and d[st.top()] <= d[j]) 
                st.pop();
            r[j] = st.empty() ? m : st.top();
            st.push(j); 
        }
        for (int j = 0; j < m; j++) 
            ans = max(ans, (i - d[j]) * (r[j] - l[j] - 1));
    }
    return ans;
}
