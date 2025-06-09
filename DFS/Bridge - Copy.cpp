#include <bits/stdc++.h>
using namespace std;

typedef long long int ll;
const int N = 1e5 + 10, INF = 1e9 + 10;

bool mark[N];
int height[N], dp[N];
vector<int> adj[N];
vector<pair<int, int>> bridge;

void find_bridge(int u, int par) {
    mark[u] = true;
    dp[u] = height[u];

    for (int v : adj[u]) {
        if (!mark[v]) {
            height[v] = height[u] + 1;
            find_bridge(v, u);
            if (dp[v] > height[u]) {
                bridge.push_back({u, v});
            }
            dp[u] = min(dp[u], dp[v]);
        }
        else if (v != par) {
            dp[u] = min(dp[u], height[v]);
        }
    }
}
