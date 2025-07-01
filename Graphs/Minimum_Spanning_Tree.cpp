#include <bits/stdc++.h>
using namespace std;

/*============================================================================================================
Description:
  Minimum Spanning Tree (MST) algorithms: Prim’s (dense & sparse variants) and Kruskal’s

  • Prim’s Algorithm
    – Version 1 (Dense graphs)
      • Uses adjacency matrix `adj[n][n]`
      • At each step, scans all n vertices to pick the next one
      • Time Complexity: O(n²)
    – Version 2 (Sparse graphs)
      • Uses adjacency list `g` of Edge{to, w}
      • Maintains a min‑priority queue (set) of candidate edges
      • Time Complexity: O(m·log(n))

  • Kruskal’s Algorithm
    – Treats the graph as a collection of edges
    – Sorts all edges by weight in O(m·log(m))
    – Uses a Disjoint Set Union (DSU) / Union–Find structure to add edges
      without forming cycles
    – Each `union` / `find` operation is essentially O(α(n)) ≈ O(1)
    – Overall Time Complexity: O(m·log(m)) ≈ O(m·log(n)) for m edges, n vertices


Properties of the minimum spanning tree:
  – A minimum spanning tree (MST) is unique if all edge weights are distinct
    Otherwise, multiple MSTs may exist, Specific algorithms usually output one MST
  – MST is also the tree with the minimum product of edge weights
  – The maximum edge weight in the MST is the minimum possible among all spanning trees
  – The maximum spanning tree can be found by negating all edge weights and applying an MST algorithm

Applications:
  – Network design (cables, roads)
  – Clustering and approximation algorithms
  – Any scenario requiring the lowest‑cost spanning structure
============================================================================================================*/


struct Prim {
    struct Edge {
        int to = -1, w = 1e9 + 10;
        bool operator<(Edge const &other) const {
            return make_pair(w, to) < make_pair(other.w, other.to);
        }
    };

    int n, INF = 1e9 + 10;
    int adj[1000][1000];
    vector<vector<Edge>> g;
    vector<pair<int, int>> MST_edges;

    bool prim() {
        int total_weight = 0;
        vector<bool> selected(n, false);
        vector<Edge> min_edge(n);
        min_edge[0].w = 0;
        for (int _ = 0; _ < n; _++) {
            int u = -1;
            for (int i = 0; i < n; i++) if (!selected[i] and (u == -1 or min_edge[i].w < min_edge[u].w)) u = i;
            if (min_edge[u].w == INF) return false;
            selected[u] = true, total_weight += min_edge[u].w;
            if (min_edge[u].to != -1) MST_edges.push_back({u, min_edge[u].to});
            for (int v = 0; v < n; v++) if (adj[u][v] < min_edge[v].w) min_edge[v] = {v, adj[u][v]};
        }
        return true;
    }
    bool prim() {
        int total_weight = 0;
        vector<bool> selected(n, false);
        vector<Edge> min_edge(n);
        min_edge[0].w = 0;
        set<Edge> q;
        q.insert({0, 0});
        for (int _ = 0; _ < n; _++) {
            if (q.empty()) return false;
            int u = q.begin()->to;
            selected[u] = true, total_weight += q.begin()->w;
            q.erase(q.begin());
            if (min_edge[u].to != -1) MST_edges.push_back({u, min_edge[u].to});
            for (Edge e : g[u]) 
                if (!selected[e.to] and e.w < min_edge[e.to].w) {
                    q.erase({e.to, min_edge[e.to].w});
                    min_edge[e.to] = {u, e.w};
                    q.insert({e.to, e.w});
                }
        }
        return true;
    }
};

struct Kruskal {
    struct Edge {
        int u, v, w;
        bool operator<(Edge const &other) {
            return w < other.w;
        }
    };

    vector<Edge> edges;
    vector<int> par, rank;
    int n, total_weight = 0;
    vector<pair<int, int>> MST_edges;

    void make_set(int u) {
        par[u] = u, rank[u] = 0;
    }
    int find_set(int u) {
        if (u == par[u]) return u;
        return par[u] = find_set(par[u]);
    }
    void union_set(int u, int v) {
        u = find_set(u), v = find_set(v);
        if (u != v) {
            if (rank[u] < rank[v]) swap(u, v);
            par[v] = u;
            if (rank[u] == rank[v]) rank[u]++;
        }
    }
    bool kruskal() {
        par.resize(n), rank.resize(n);
        for (int u = 0; u < n; u++) make_set(u);
        sort(edges.begin(), edges.end());
        for (Edge e : edges) 
            if (find_set(e.u) != find_set(e.v)) {
                total_weight += e.w;
                MST_edges.push_back({e.u, e.v});
                union_set(e.u, e.v);
                if (MST_edges.size() == n - 1) return true;
            }
        return false;    
    }
};

/*============================================================================================================
Description:
  This code computes the minimum spanning tree (MST) of a connected undirected weighted graph,
  classifies every edge into one of:
    - appears in all MSTs 
    - appears in at least one MST 
    - appears in no MST
  and also computes the second-best MST

Key Steps:
  1. Kruskal Sweep with Bucketed Weights + Bridge Detection
    • Sort edges by weight, and process in equal-weight buckets
    • Use DSU of components formed by edges of lesser weight
    • Build an auxiliary component-level graph for the current bucket
    • Run Tarjan’s DFS to detect bridges, which mark edges that belong to all MSTs
    • Union all bucket edges into the DSU and record which edges enter the MST

  2. LCA Preprocessing (Binary Lifting) for 'max-edge-on-path' Queries
    • Root the MST at node 1 and run DFS to record depths and immediate parents
    • Build parent[k][v] and maxWeight[k][v] tables to support O(log n) queries for the maximum edge on the path to any ancestor

  3. Second-Best Candidate & 'Could-Be-In-MST' Detection
    • For each non-MST edge e(u,v):
      – Find max-weight edge f on the path between u and v in the MST
      – Compute candidate cost: `cand = MST_total + w(e) - w(f)`
      – If `cand == MST_total`, then swapping e with f yields another MST ⇒ e is in some MST
      – If `cand > MST_total`, it’s a valid second-best spanning tree candidate

Time Complexity:
  • Sorting + DSU: O(m.log(m) + m.α(n))
  • Bridge detection (over buckets): O(m + n)
  • LCA preprocessing: O(n.log(n))
  • Edge-by-edge analysis: O(m.log(n))
  • Overall: O(m.log(m))
============================================================================================================*/

struct SecondBestMST {
    struct Edge {
        int u, v, w, id;
    };
    
    struct DSU {
        vector<int> par, rank;
        
        DSU(int n): par(n + 1), rank(n + 1, 0) {
            iota(par.begin(), par.end(), 0);
        }
        int find(int u) {
            return par[u] == u ? u : par[u] = find(par[u]);
        }
        bool unite(int u, int v) {
            u = find(u), v = find(v);
            if (u == v) return false;
            if (rank[u] < rank[v]) swap(u, v);
            par[v] = u;
            if (rank[u] == rank[v]) rank[u]++;
            return true;
        }
    };
    
    vector<int> depth;
    vector<Edge> edges;
    vector<vector<tuple<int, int, int>>> adj;
    vector<vector<int>> parent, maxEdge, maxWeight;
    vector<bool> in_mst, could_be_in_MST, always_in_MST;
    int n, m, LOG, mst_weight = 0, second_weight = 1e9 + 10;
    
    void build_mst_and_classify(DSU &dsu) {
        sort(edges.begin(), edges.end(), [] (auto &a, auto &b) {return a.w < b.w;});
        int i = 0;
        while (i < m) {
            int j = i, w = edges[i].w, k = 0;
            while (j < m and edges[j].w == w) j++;
            vector<Edge> bucket(edges.begin() + i, edges.begin() + j);
            unordered_map<int, int> comp_id;
            for (Edge &e : bucket) {
                int cu = dsu.find(e.u), cv = dsu.find(e.v);
                if (cv == cu) continue;
                if (!comp_id.count(cu)) comp_id[cu] = k++;
                if (!comp_id.count(cv)) comp_id[cv] = k++;
            }
            vector<vector<pair<int, int>>> H(k);
            for (Edge &e : bucket) {
                int cu = dsu.find(e.u), cv = dsu.find(e.v);
                if (cv == cu) continue;
                int iu = comp_id[cu], iv = comp_id[cv];
                H[iu].emplace_back(iv, e.id);
                H[iv].emplace_back(iu, e.id);
            }
            vector<int> disc(k, -1), low(k);
            function <void (int, int)> DFS = [&] (int u, int peid) {
                static int t;
                if (peid < 0) t = 0;
                disc[u] = low[u] = t++;
                for (auto &pr : H[u]) {
                    int v = pr.first, eid = pr.second;
                    if (eid == peid) continue;
                    if (disc[v] < 0) {
                        DFS(v, eid);
                        low[u] = min(low[u], low[v]);
                        if (low[v] > disc[u]) always_in_MST[eid] = true;
                    }
                    else low[u] = min(low[u], disc[v]);
                }
            };
            for (int u = 0; u < k; u++) if (disc[u] < 0) DFS(u, -1);
            for (Edge &e : bucket) {
                if (dsu.unite(e.u, e.v)) {
                    in_mst[e.id] = true;
                    mst_weight += e.w;
                    adj[e.u].emplace_back(e.v, e.w, e.id);
                    adj[e.v].emplace_back(e.u, e.w, e.id);
                }
            }
            i = j;
        }
    }
    void preprocess_lca() {
        LOG = 1; while ((1 << LOG) <= n) LOG++;
        depth.assign(n + 1, 0);
        parent.assign(LOG, vector<int>(n + 1, 0));
        maxWeight.assign(LOG, vector<int>(n + 1, 0));
        maxEdge.assign(LOG, vector<int>(n + 1, -1));
        function <void (int, int)> DFS = [&] (int u, int p) {
            for (auto [v, w, id] : adj[u]) {
                if (v == p) continue;
                depth[v] = depth[u] + 1, parent[0][v] = u, maxWeight[0][v] = w, maxEdge[0][v] = id;
                DFS(v, u);
            }
        };
        DFS(1, 0);
        for (int k = 1; k < LOG; k++)
            for (int u = 1; u <= n; u++) {
                int mid = parent[k - 1][u];
                parent[k][u] = parent[k - 1][mid];
                if (maxWeight[k - 1][u] >= maxWeight[k - 1][mid]) maxWeight[k][u] = maxWeight[k - 1][u], maxEdge[k][u] = maxEdge[k - 1][u];
                else maxWeight[k][u] = maxWeight[k - 1][mid], maxEdge[k][u] = maxEdge[k - 1][mid];
            }
    }
    pair<int, int> get_max_on_path(int u, int v) {
        if (depth[u] < depth[v]) swap(u, v);
        int best_w = 0, best_id = -1;
        for (int k = LOG - 1; k >= 0; k--)
            if (depth[u] - (1 << k) >= depth[v]) {
                if (maxWeight[k][u] > best_w) best_w = maxWeight[k][u], best_id = maxEdge[k][u];
                u = parent[k][u];
            }
        if (u == v) return {best_w, best_id};
        for (int k = LOG - 1; k >= 0; k--)
            if (parent[k][u] != parent[k][v]) {
                if (maxWeight[k][u] > best_w) best_w = maxWeight[k][u], best_id = maxEdge[k][u];
                if (maxWeight[k][v] > best_w) best_w = maxWeight[k][v], best_id = maxEdge[k][v];
                u = parent[k][u], v = parent[k][v];
            }
        if (maxWeight[0][u] > best_w) best_w = maxWeight[0][u], best_id = maxEdge[0][u];
        if (maxWeight[0][v] > best_w) best_w = maxWeight[0][v], best_id = maxEdge[0][v];
        return {best_w, best_id};
    }
    void compute_second() {
        for (Edge &e : edges) {
            if (in_mst[e.id]) {could_be_in_MST[e.id] = true; continue;}
            auto pr = get_max_on_path(e.u, e.v);
            int cand = mst_weight + e.w - pr.first;
            if (cand == mst_weight) could_be_in_MST[e.id] = true;
            else if (cand > mst_weight) second_weight = min(second_weight, cand);
        }
    }
    void init() {
        cin >> n >> m;
        edges.resize(m);
        for (int i = 0; i < m; i++) {
            cin >> edges[i].u >> edges[i].v >> edges[i].w;
            edges[i].id = i;
        }
        in_mst.assign(m, false), could_be_in_MST.assign(m, false), always_in_MST.assign(m, false), adj.assign(n + 1, {});
        DSU dsu = DSU(n);
        build_mst_and_classify(dsu);
        preprocess_lca();
        compute_second();
    }
};

/*============================================================================================================
Description:
  Count the number of spanning trees in a connected undirected graph (with multiple edges/loops)
  using Kirchhoff’s Matrix‑Tree Theorem

Definitions:
  • A (adjacency matrix):
    A[u][v] = number of edges between u and v
    Note: a loop at u contributes twice to A[u][u]
  • D (degree matrix):
    Diagonal matrix where D[u][u] = degree of vertex u = ∑ₖ A[u][k]
  • Laplacian L = D – A

  • Build L in O(n²)
  • Form any cofactor by removing one row & one column → an (n–1)×(n–1) minor
  • Compute its determinant via Gaussian elimination in O(n³)
  • The absolute value of that determinant = number of spanning trees

Overall Time Complexity: O(n³)
============================================================================================================*/