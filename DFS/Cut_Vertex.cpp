#include <bits/stdc++.h>
using namespace std;

typedef long long int ll;
const int N = 1e5 + 10, INF = 1e9 + 10;

int root, height[N], dp[N];
bool mark[N], is_cut[N];
vector<int> adj[N];

void find_cut(int u, int par) {
    mark[u] = true;
    dp[u] = height[u];
    int child_cnt = 0;
    for (int v : adj[u]) {
        if (!mark[v]) {
            height[v] = height[u] + 1;
            find_cut(v, u);
            if (u != root && dp[v] >= height[u]) {
                is_cut[u] = true;
            }
            dp[u] = min(dp[u], dp[v]);
            child_cnt++;
        }
        else if (v != par) {
            dp[u] = min(dp[u], height[v]);
        }
    }
    if (u == root && child_cnt > 1) {
        is_cut[u] = true;
    }
}
