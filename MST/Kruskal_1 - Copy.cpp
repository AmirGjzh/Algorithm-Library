#include <bits/stdc++.h>
using namespace std;

typedef long long int ll;
const int N = 1e5 + 10, INF = 1e9 + 10;

int n, m, par[N];
vector<pair<int, pair<int, int>>> edges;

int root(int u) {
    if (par[u] == u) {
        return u;
    }
    return par[u] = root(par[u]);
}

void Union(int u, int v) {
    par[root(u)] = root(v);
}

void kruskal() {
    sort(edges.begin(), edges.end());

    for (int u = 0; u < n; u++) {
        par[u] = u;
    }

    ll MST_LENGTH = 0;

    for (int e = 0; e < m; e++) {
        int u = edges[e].second.first;
        int v = edges[e].second.second;

        if (root(u) != root(v)) {
            MST_LENGTH += edges[e].first;
            Union(u, v);
        }
    }
}
