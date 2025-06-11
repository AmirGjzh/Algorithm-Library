#include <bits/stdc++.h>
using namespace std;

typedef long long int ll;
const int MOD = 1e9 + 7;
const int N = 1e4 + 10, INF = 1e9 + 10;

bool mark[N];
int match_size, match[N * 2];
set<int> up, down;
vector<int> adj[N];

bool dfs(int u) {
    mark[u] = true;
    for (int v : adj[u]) {
        if (match[v] == -1 || (!mark[match[v]] && dfs(match[v]))) {
            match[u] = v;
            match[v] = u;
            return true;
        }
    }
    return false;
}

void find_matching() {
    fill(match, match + N * 2, -1);

    for (int u : up) {
        fill(mark, mark + N, false);
        if (dfs(u)) {
            match_size++;
        }
    }
}
