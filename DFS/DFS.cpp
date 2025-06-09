#include <bits/stdc++.h>
using namespace std;

typedef long long int ll;
const int N = 1e5 + 10, INF = 1e9 + 10;

bool mark[N];
vector<int> adj[N];
int time_now, color[N], start[N], finish[N];

void dfs(int u) {
    mark[u] = true;
    color[u] = 1;
    start[u] = time_now++;
    for (int v : adj[u]) {
        if (!mark[v]) {
            dfs(v);
        }
    }
    color[u] = 2;
    finish[u] = time_now;
}
