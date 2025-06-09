#include <bits/stdc++.h>
using namespace std;

typedef long long int ll;
const int N = 1e5 + 10, INF = 1e9 + 10;

int dis[N], par[N];
bool in_queu[N];
vector<pair<int, int>> adj[N];

void spfa(int s, int e) {
    fill(dis, dis + N, INF);
    fill(par, par + N, -1);

    queue<int> q;
    dis[s] = 0;
    in_queu[s] = true;
    q.push(s);

    while (!q.empty()) {
        int u = q.front();
        in_queu[u] = false;
        q.pop();

        if (dis[u] > dis[e]) {
            continue;
        }

        for (auto v : adj[u]) {
            if (dis[u] + v.second < dis[v.first]) {
                dis[v.first] = dis[u] + v.second;
                par[v.first] = u;
                if (!in_queu[v.first]) {
                    q.push(v.first);
                    in_queu[v.first] = true;
                }
            }
        }
    }
}
