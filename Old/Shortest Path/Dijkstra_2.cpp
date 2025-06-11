#include <bits/stdc++.h>
using namespace std;

typedef long long int ll;
const int N = 1e5 + 10, INF = 1e9 + 10;

int n, dis[N], par[N], adj[N][N];
bool mark[N];

void dijkstra(int s) {
    fill(dis, dis + N, INF);
    fill(par, par + N, -1);
    dis[s] = 0;

    for (int i = 0; i < n; i++) {
        int u = -1;
        int du = INF;

        for (int v = 0; v < n; v++) {
            if (!mark[v] && dis[v] <= du) {
                du = dis[v];
                u = v;
            }
        }
        mark[u] = true;
        for (int v = 0; v < n; v++) {
            if (adj[u][v] != INF && !mark[v]) {
                if (dis[v] > dis[u] + adj[u][v]) {
                    dis[v] = dis[u] + adj[u][v];
                    par[v] = u;
                }
            }
        }
    }
}
