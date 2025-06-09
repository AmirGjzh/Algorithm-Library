#include <bits/stdc++.h>
using namespace std;

typedef long long int ll;
const int N = 1e5 + 10, INF = 1e9 + 10;

int dis[N];
bool mark[N];
vector<pair<int, int>> adj[N];

void bfs(int s) {
    fill(dis, dis + N, INF);
    dis[s] = 0;
    deque<int> q;
    q.push_back(s);
    while (!q.empty()) {
        int u = q.front();
        q.pop_front();
        if (!mark[u]) {
            mark[u] = true;
            for (auto v : adj[u]) {
                if (dis[u] + v.second < dis[v.first]) {
                    dis[v.first] = dis[u] + v.second;
                    if (v.second == 0) {
                        q.push_front(v.first);
                    }
                    else {
                        q.push_back(v.first);
                    }
                }
            }
        }
    }
}
