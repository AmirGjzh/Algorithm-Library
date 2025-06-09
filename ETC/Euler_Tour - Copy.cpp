#include <bits/stdc++.h>
using namespace std;

typedef long long int ll;
const int N = 1e5 + 10, M = 1e6 + 10, INF = 1e9 + 10;

int pt[N];
bool mark[M];
vector<int> ans;
vector<pair<int, int>> adj[N];

void euler_tour(int u) {
    while (pt[u] < adj[u].size()) {
        int v = adj[u][pt[u]].first;
        int edge = adj[u][pt[u]].second;
        pt[u]++;

        if (!mark[edge]) {
            mark[edge] = true;
            euler_tour(v);
        }
    }
    ans.push_back(u);
}
