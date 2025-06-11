#include <bits/stdc++.h>
using namespace std;

typedef long long int ll;
const int N = 1e5 + 10, INF = 1e9 + 10;

int n, dis[N], par[N];
vector<pair<int, int>> adj[N];

void dijkstra(int s) {
    fill(dis, dis + n, INF);
    fill(par, par + n, -1);
    
    dis[s] = 0;
    set<pair<int, int>> q;
    for (int i = 0; i < n; i++) {
        q.insert({dis[i], i});
    }

    while (!q.empty()) {
        int u = q.begin() -> second, d = q.begin() -> first;
        q.erase(q.begin());
        for (auto v : adj[u]) {
            if (dis[v.first] > dis[u] + v.second) {
                q.erase({dis[v.first], v.first});
                dis[v.first] = dis[u] + v.second;
                par[v.first] = u;
                q.insert({dis[v.first], v.first});
            }
        }
    }
}
