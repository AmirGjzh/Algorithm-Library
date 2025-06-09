#include <bits/stdc++.h>
using namespace std;

typedef long long int ll;
const int N = 1e5 + 10, INF = 1e9 + 10;

int n, dis[N][N];

void floyd_warshall() {
    for (int k = 0; k < n; k++) {
        for (int u = 0; u < n; u++) {
            for (int v = 0; v < n; v++) {
                dis[u][v] = min(dis[u][v], dis[u][k] + dis[k][v]);
            }
        }
    }
}
