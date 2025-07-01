#include <bits/stdc++.h>
using namespace std;

/*============================================================================================================
Dijkstra (dense, set-based, and heap-based versions)

Description:
  • Finds shortest paths from a source in a graph with non-negative weights
  • Three implementations: O(n²) naive, O(m.log(n)) via set, and O(m.log(n)) via min-heap

Applications:
  • Routing in positive-weighted networks
  • Subroutine for Johnson's algorithm & other combinatorial tasks
  • Basic heuristic for A* (with zero heuristic)

Properties:
  • When a node u is popped from the priority queue for the i-th time, the distance d at that moment
    is the i-th smallest shortest-path distance to u, This underpins the k-th shortest paths algorithm

Notes:
  • Priority_queue vs set performance:
    Though both `priority_queue` and `set` implementations have O(m.log(n)) theoretical complexity, 
    `priority_queue` is significantly faster in practice

Order:
  • Dense: O(n²)
  • Set & heap: O(m.log(n))
============================================================================================================*/

struct Dijkstra {
    int n, INF = 1e9 + 10;
    vector<int> dis, par;
    vector<vector<pair<int, int>>> g;

    void dijkstra_basic(int s) {
        dis.assign(n, INF), par.assign(n, -1);
        vector<bool> used(n, false);
        dis[s] = 0;
        for (int _ = 0; _ < n; _++) {
            int u = -1;
            for (int i = 0; i < n; i++) if (!used[i] and (u == -1 or dis[i] < dis[u])) u = i;
            if (dis[u] == INF) break;
            used[u] = true;
            for (auto [v, w] : g[u]) if (dis[v] > dis[u] + w) {dis[v] = dis[u] + w, par[v] = u;}
        }
    }
    void dijkstra_set(int s) {
        dis.assign(n, INF), par.assign(n, -1);
        dis[s] = 0;
        set<pair<int, int>> q;
        q.insert({dis[s], s});
        while (q.size()) {
            int u = q.begin()->second;
            q.erase(q.begin());
            for (auto [v, w] : g[u]) 
                if (dis[v] > dis[u] + w) {
                    q.erase({dis[v], v});
                    dis[v] = dis[u] + w;
                    par[v] = u;
                    q.insert({dis[v], v});
                }
        }
    }
    void dijkstra_pq(int s) {
        dis.assign(n, INF), par.assign(n, -1);
        priority_queue<pair<int, int>, vector<pair<int, int>>, greater<pair<int, int>>> q;
        dis[s] = 0;
        q.push({dis[s], s});
        while (q.size()) {
            int u = q.top().second, d = q.top().first;
            q.pop();
            if (d != dis[u]) continue;
            for (auto [v, w] : g[u]) 
                if (dis[v] > dis[u] + w) {
                    dis[v] = dis[u] + w;
                    par[v] = u;
                    q.push({dis[v], v});
                }
        }
    }
    vector<int> restore_path(int t) {
        vector<int> path;
        for (int u = t; u != -1; u = par[u]) path.push_back(u);
        reverse(path.begin(), path.end());
        return path;
    }
};

/*============================================================================================================
Bellman-Ford Algorithm

Description:
  • Finds single-source shortest paths in graphs with negative weights (no negative cycle)
  • Detects negative cycles reachable from source

Applications:
  • Currency arbitrage, graphs with mixed weights
  • Preprocessing for Johnson’s reweighting

Notes:
  • Runs n iterations over all edges. Early exit if no relaxation occurs
  • Negative cycle detection: if relaxation still possible on nth iteration → cycle present
  • You capture and print one cycle via parent-pointer tracing
  • To find any negative cycle (not necessarily reachable from a given source), you can
    initialize all distances to zero and run Bellman–Ford from “every” vertex at once
    (i.e set `dis.assign(n, 0)` before relaxation)
  • In a directed graph, to detect a negative cycle on some path from 's' to 't', after the final Bellman–Ford relaxation, 
    check whether any vertex relaxed in that last pass can reach 't'

Order: O(n·m)
============================================================================================================*/

struct BellmanFord {
    struct Edge {
        int u, v, w;
    };

    vector<Edge> edges;
    vector<int> dis, par;
    int n, m, INF = 1e9 + 10;

    void bellman_ford(int s) {
        dis.assign(n, INF), par.assign(n, -1);
        dis[s] = 0;
        int x;
        for (int _ = 0; _ < n; _++) { 
            x = -1;
            for (Edge e : edges)
                if (dis[e.u] < INF and dis[e.v] > dis[e.u] + e.w) {
                    dis[e.v] = max(-INF, dis[e.u] + e.w);
                    par[e.v] = e.u;
                    x = e.v;
                }    
            if (x == -1) break;          
        }
        if (x != -1) {
            for (int i = 0; i < n; i++) x = par[x];
            vector<int> cycle;
            for (int cur = x; ; cur = par[cur]) {
                cycle.push_back(cur);
                if (cur == x and cycle.size() > 1) break;
            }
            reverse(cycle.begin(), cycle.end());
        }
    }
    vector<int> restore_path(int t) {
        vector<int> path;
        for (int u = t; u != -1; u = par[u]) path.push_back(u);
        reverse(path.begin(), path.end());
        return path;
    }
};

/*============================================================================================================
Shortest Path Faster Algorithm (SPFA)

Description:
  • Optimized Bellman-Ford using a queue, effective in sparse graphs

Applications:
  • Mixed-weight graphs where Dijkstra can't be used

Notes:
  • Tracks `cnt[v]` (number of times v was enqueued). If `cnt[v] > n`, a negative cycle exists
  • May perform slowly or loop on adversarial graphs—even exponential in worst cases
  
Order: O(n·m) worst, much faster in practice
============================================================================================================*/

struct SPFA {
    vector<int> dis, par;
    int n, INF = 1e9 + 10;
    vector<vector<pair<int, int>>> g;

    bool spfa(int s) {
        dis.assign(n, INF), par.assign(n, -1);
        vector<int> cnt(n, 0);
        vector<bool> inqueue(n, false);
        queue<int> q;
        dis[s] = 0;
        q.push(s);
        inqueue[s] = true;
        while (q.size()) {
            int u = q.front();
            q.pop();
            inqueue[u] = false;
            for (auto [v, w] : g[u]) 
                if (dis[v] > dis[u] + w) {
                    dis[v] = dis[u] + w, par[v] = u;
                    if (!inqueue[v]) {
                        q.push(v);
                        inqueue[v] = true, cnt[v]++;
                        if (cnt[v] > n) return false;
                    }
                }
        }
        return true;
    }
    vector<int> restore_path(int t) {
        vector<int> path;
        for (int u = t; u != -1; u = par[u]) path.push_back(u);
        reverse(path.begin(), path.end());
        return path;
    }
};

/*============================================================================================================
0–1 BFS & Dial’s Algorithm

Description:
  • 0–1 BFS: Finds shortest paths in graphs where edge weights are either 0 or 1
  • Dial’s algorithm: Extends this idea to integer edge weights in the range [0…k], using bucketed queues

Applications:
  1. Binary-cost routing problems (e.g, 0–1 BFS: lab maze, weighted grids)
  2. Dial’s algorithm: fast shortest paths when edge weights are bounded small integers
  3. Dial’s algorithm is often used in networking and discrete event simulations (bucket queue)

Notes:
  • 0–1 BFS: Use a deque, push_front for 0-weight edges, push_back for 1-weight. Runs in O(n + m)
  • Dial’s algorithm:
    – If all edge weights ≤ k, maintain k+1 buckets (arrays or lists), each holding vertices grouped by tentative distance modulo (k+1)
    – Process buckets cyclically: when the current bucket is empty, move to the next non-empty—ensuring you always process the smallest unreached distance
    – Time complexity: O(m + n·k), which beats O(m.log(n)) for small k
  • Dial’s method is sometimes called a "bucket queue" or "bounded-height priority queue"

Order:
  • 0–1 BFS: O(n + m)  
  • Dial’s algorithm: O(m + n·k)
============================================================================================================*/

struct BFS_01 {
    vector<int> dis;
    int n, INF = 1e9 + 10;
    vector<vector<pair<int, int>>> g;

    void bfs_01(int s) {
        dis.assign(n, INF);
        dis[s] = 0;
        deque<int> q;
        q.push_front(s);
        while (q.size()) {
            int u = q.front();
            q.pop_front();
            for (auto [v, w] : g[u]) 
                if (dis[v] > dis[u] + w) {
                    dis[v] = dis[u] + w;
                    if (w == 1) q.push_back(v);
                    else q.push_front(v);
                }
        }
    }
};

/*============================================================================================================
D’Esopo–Pape Algorithm

Description:
  • Like SPFA but deque-based with push-front/back and state `m[v]`
  • Faster in many practical scenarios for negative or mixed graphs (no cycles)

Applications:
  • General shortest-path with possible negative edges, when SPFA is too slow

Notes:
  • Holds “status” array `m[] = {0,1,2}` to manage deque insertion
  • Worst-case time is exponential—avoid on worst-case tests

Order: Exponential worst-case, often much faster
============================================================================================================*/

struct DEsopo_Pape {
    vector<int> dis, par;
    int n, INF = 1e9 + 10;
    vector<vector<pair<int, int>>> g;

    void solve(int s) {
        dis.assign(n, INF), par.assign(n, -1);
        vector<int> m(n, 2);
        deque<int> q;
        dis[s] = 0;
        q.push_back(s);
        while (q.size()) {
            int u = q.front();
            q.pop_front();
            m[u] = 0;
            for (auto [v, w] : g[u]) 
                if (dis[v] > dis[u] + w) {
                    dis[v] = dis[u] + w, par[v] = u;
                    if (m[v] == 2) {m[v] = 1; q.push_back(v);}
                    else if (m[v] == 0) {m[v] = 1; q.push_front(v);}
                }
        }
    }
    vector<int> restore_path(int t) {
        vector<int> path;
        for (int u = t; u != -1; u = par[u]) path.push_back(u);
        reverse(path.begin(), path.end());
        return path;
    }
};

/*============================================================================================================
Floyd–Warshall (All-Pairs Shortest Paths)

Description:
  • Computes shortest paths between all pairs in a weighted directed graph
  • Supports negative weights, but not negative cycles

Applications:
  • Dense graphs (n ≤ ~500)
  • Known uses: APSP, path-counting, graph centrality, triangle detection

Notes:
  • Initialize `dis[i][i] = 0`, `dis[i][j] = weight` or INF
  • Handles loops correctly if initialized

Order: O(n³)
============================================================================================================*/

/*============================================================================================================
Floyd–Warshall (All-Pairs Shortest Paths)

Description:
  • Computes shortest paths between all pairs in a weighted directed graph
  • Supports negative weights, but not negative cycles

Applications:
  • Dense graphs (n ≤ ~500)
  • Known uses: APSP, path-counting, graph centrality, triangle detection

Notes:
  • Initialize `dis[i][i] = 0`, `dis[i][j] = weight` or INF
  • Handles loops correctly if initialized
  • `find_negative_cycle()` marks any pair (u, v) as -INF if there exists a negative cycle 
    reachable on some path from u to v, indicating no well-defined shortest path

Order: O(n³)
============================================================================================================*/


struct FloydWarshall {
    int dis[500][500];
    int n, INF = 1e9 + 10;

    void floyd_warshall() {
        for (int k = 0; k < n; k++)
            for (int u = 0; u < n; u++) 
                for (int v = 0; v < n; v++)
                    if (dis[u][k] < INF and dis[k][v] < INF)
                        dis[u][v] = min(dis[u][v], dis[u][k] + dis[k][v]);
    }
    void find_negative_cycle() {
        floyd_warshall();
        for (int u = 0; u < n; u++)
            for (int v = 0; v < n; v++) 
                for (int t = 0; t < n; t++)
                    if (dis[u][t] < INF and dis[t][v] < INF and dis[t][t] < 0) 
                        dis[u][v] = -INF;
    }
};

/*============================================================================================================
k‑Length Path Counts and Shortest Paths via Matrix Exponentiation

1. Count Paths of fixed length:
   • Raise adjacency matrix G to the power k; entry (i,j) gives # paths of length k
   • Works even with loops or multi-edges (counts via multiplicity in G)

2. Shortest path with exactly k edges:
   • Use min-plus matrix multiplication: C[i][j] = min over p (A[i][p] + B[p][j])
   • Perform binary exponentiation for O(n³·log(k))

3. Up to k edges:
   • Duplicate nodes v → v′, add zero-weight edges v → v′ and v′ → v′
   • Now paths of length ≤ k in original equal paths of length exactly k+1 ending at v′

Order: O(n³·log(k))
============================================================================================================*/
