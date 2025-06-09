#include <bits/stdc++.h>
using namespace std;

typedef long long int ll;
const int N = 1e5 + 10, INF = 1e9 + 10;

int n, c, col[N];
bool mark[N];
stack<int> sources; 
vector<int> adj[N], rev_adj[N];

void scc(int u) {
    mark[u] = true;
    col[u] = c;
    for (int v : rev_adj[u]) {
        if (!mark[v]) {
            scc(v);
        }
    }
}

void dfs(int u) {
    mark[u] = true;
    for (int v : adj[u]) {
        if (!mark[v]) {
            dfs(v);
        }
    }
    sources.push(u);
}

void find_scc() {
    for (int u = 0; u < n; u++) {
        if (!mark[u]) {
            dfs(u);
        }
    }
    fill(mark, mark + N, false);
    while (!sources.empty()) {
        int u = sources.top();
        sources.pop();
        if (!mark[u]) {
            scc(u);
            c++;
        }
    }
}
