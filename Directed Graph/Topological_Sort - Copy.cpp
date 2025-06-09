#include <bits/stdc++.h>
using namespace std;

typedef long long int ll;
const int N = 1e5 + 10, INF = 1e9 + 10;

int n, deg[N];
set<int> sources;
vector<int> ans, adj[N];

bool topological_sort() {
    for (int u = 0; u < n; u++) {
        if (deg[u] == 0) {
            sources.insert(u);
        }
    }

    for (int i = 0; i < n; i++) {
        if (sources.empty()) {
            return false;
        }

        int u = *sources.begin();
        sources.erase(sources.begin());
        
        for (int v : adj[u]) {
            if (--deg[v] == 0) {
                sources.insert(v);
            }
        }
        ans.push_back(u);
    }
    
    return true;
}
