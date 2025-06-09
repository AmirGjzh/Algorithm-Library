#include <bits/stdc++.h>
using namespace std;

typedef long long int ll;
const int N = 1e5 + 10, INF = 1e9 + 10;

int n, dis[N], par[N];
vector<pair<int, int>> edge;
vector<int> weight;

bool bellman_ford(int s) {
    fill(dis, dis + N, INF);
    fill(par, par + N, -1);

    dis[s] = 0;
    for (int i = 0; i < n - 1; i++) {
        for (int e = 0; e < edge.size(); e++) {
            int u = edge[e].first;
            int v = edge[e].second;

            if (dis[v] > dis[u] + weight[e]) {
                dis[v] = dis[u] + weight[e];
                par[v] = u;
            }
        }
    }

    for (int e = 0; e < edge.size(); e++) {
        int u = edge[e].first;
        int v = edge[e].second;

        if (dis[v] > dis[u] + weight[e]) {
            return false;
        }
    }
}
