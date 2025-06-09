#include <bits/stdc++.h>
using namespace std;

typedef long long int ll;
const int N = 1e5 + 10, INF = 1e9 + 10;

int n, m, tree_id[N];
vector<pair<int, pair<int, int>>> edges;

void kruskal() {
    sort(edges.begin(), edges.end());

    for (int u = 0; u < n; u++) {
        tree_id[u] = u;
    }

    ll MST_LENGTH = 0;

    for (int e = 0; e < m; e++) {
        int u = edges[e].second.first;
        int v = edges[e].second.second;

        if (tree_id[u] != tree_id[v]) {
            MST_LENGTH += edges[e].first;

            int old_id = tree_id[v];
            for (int r = 0; r < n; r++) {
                if (tree_id[r] == old_id) {
                    tree_id[r] = tree_id[u];
                }
            }
        }
    }
}
