#include <bits/stdc++.h>
using namespace std;

typedef long long int ll;
const int N = 1e5 + 10, INF = 1e9 + 10;

int n, min_edge[N];
bool mark[N];
vector<pair<int, int>> adj[N];

void prim(int u) {
    fill(min_edge, min_edge + N, INF);

    min_edge[u] = 0;
    ll MST_LENGTH = 0;

    for (int i = 0; i < n; i++) {
        int v = -1;
        for (int r = 0; r < n; r++) {
            if (!mark[r] && (v == -1 || min_edge[r] < min_edge[v])) {
                v = r;
            }
        }
        
        mark[v] = true;
        MST_LENGTH += min_edge[v];

        for (auto r : adj[v]) {
            if (!mark[r.first]) {
                min_edge[r.first] = min(min_edge[r.first], r.second);
            }    
        }
    }
}
