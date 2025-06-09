#include <bits/stdc++.h>
using namespace std;

typedef long long int ll;
const int N = 1e5 + 10, INF = 1e9 + 10;

int n, c;
bool ans[N];
map<pair<int, bool>, int> col;
map<pair<int, bool>, bool> mark;
stack<pair<int, bool>> sources; 
map<pair<int, bool>, vector<pair<int, bool>>> adj, rev_adj;

void scc(pair<int, bool> u) {
    mark[u] = true;
    col[u] = c;
    for (auto v : rev_adj[u]) {
        if (!mark[v]) {
            scc(v);
        }
    }
}

void dfs(pair<int, bool> u) {
    mark[u] = true;
    for (auto v : adj[u]) {
        if (!mark[v]) {
            dfs(v);
        }
    }
    sources.push(u);
}

void find_scc() {
    for (int u = 0; u < n; u++) {
        if (!mark[{u, true}]) {
            dfs({u, true});
        }
        if (!mark[{u, false}]) {
            dfs({u, false});
        }
    }
    mark.clear();
    while (!sources.empty()) {
        auto u = sources.top();
        sources.pop();
        if (!mark[u]) {
            scc(u);
            c++;
        }
    }
}

bool two_sat() {
    find_scc();
    for (int u = 0; u < n; u++) {
        ans[u] = col[{u, true}] > col[{u, false}];
        if (col[{u, true}] == col[{u, false}]) {
            return false;
        }
    }
    return true;
}
