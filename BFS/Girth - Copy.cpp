#include <bits/stdc++.h>
using namespace std;

typedef long long int ll;
const int N = 1e5 + 10, INF = 1e9 + 10;

int n, dis[N], par[N];
vector<int> adj[N];

int bfs(int s, int max_depth) {
    fill(dis, dis + N, INF);
    dis[s] = 0;
    par[s] = -1;
    queue<int> q;
    q.push(s);
    while (!q.empty()) {
        int u = q.front();
        q.pop();
        if (dis[u] == max_depth) {
            return INF;
        }
        for (int v : adj[u]) {
            if (dis[v] > dis[u] + 1) {
                dis[v] = dis[u] + 1;
                par[v] = u;
                q.push(v);
            }
            else if (v != par[u]) {
                return dis[v] + dis[u] + 1;
            }
        }
    }
    return INF;
}

int find_min_cycle() {
    int min_cycle = INF;
    for (int u = 0; u < n; u++) {
        min_cycle = min(min_cycle, bfs(u, min_cycle / 2));
    }
    return min_cycle;
}
