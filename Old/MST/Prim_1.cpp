#include <bits/stdc++.h>
using namespace std;

typedef long long int ll;
const int N = 1e5 + 10, INF = 1e9 + 10;

int min_edge[N];
bool mark[N];
vector<pair<int, int>> adj[N];

void prim(int u) {
    fill(min_edge, min_edge + N, INF);

    set<pair<int, int>> s;
    s.insert({0, u});
    min_edge[u] = 0;

    ll MST_LENGTH = 0;

    while (!s.empty()) {
        int v = s.begin() -> second;
        s.erase(s.begin());
        mark[v] = true;

        MST_LENGTH += min_edge[v];

        for (auto r : adj[v]) {
            if (!mark[r.first] && min_edge[r.first] > r.second) {
                s.erase({min_edge[r.first], r.first});
                min_edge[r.first] = r.second;
                s.insert({min_edge[r.first], r.first});
            }
        }
    }
}
