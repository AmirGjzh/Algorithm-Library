#include <bits/stdc++.h>
using namespace std;
using ld = long double;
using ll = long long int;

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
        int to = -1; ll w = LLONG_MAX;
        bool operator<(Edge const &other) const {
            return make_pair(w, to) < make_pair(other.w, other.to);
        }
    };

    int n;
    ll adj[1000][1000];
    vector<vector<Edge>> G;
    vector<pair<int, int>> MST_edges;

    bool prim() {
        ll total_weight = 0;
        vector<bool> selected(n, false);
        vector<Edge> min_edge(n); min_edge[0].w = 0;
        for (int _ = 0; _ < n; _++) {
            int u = -1;
            for (int i = 0; i < n; i++) if (!selected[i] and (u == -1 or min_edge[i].w < min_edge[u].w)) u = i;
            if (min_edge[u].w == LLONG_MAX) return false;
            selected[u] = true, total_weight += min_edge[u].w;
            if (min_edge[u].to != -1) MST_edges.push_back({u, min_edge[u].to});
            for (int v = 0; v < n; v++) if (adj[u][v] < min_edge[v].w) min_edge[v] = {v, adj[u][v]};
        }
        return true;
    }
    bool prim() {
        ll total_weight = 0;
        vector<bool> selected(n, false);
        vector<Edge> min_edge(n);
        min_edge[0].w = 0; set<Edge> q; q.insert({0, 0});
        for (int _ = 0; _ < n; _++) {
            if (q.empty()) return false;
            int u = q.begin()->to;
            selected[u] = true, total_weight += q.begin()->w;
            q.erase(q.begin());
            if (min_edge[u].to != -1) MST_edges.push_back({u, min_edge[u].to});
            for (Edge e : G[u]) 
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
        int u, v; ll w;
        bool operator<(Edge const &other) {return w < other.w;}
    };

    vector<Edge> edges;
    vector<int> par, rank;
    int n; ll total_weight = 0;
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
        for (Edge &e : edges) 
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
This `MST` struct computes the Minimum Spanning Tree (MST) of a connected undirected weighted graph and
provides rich edge classification and replacement information, including:
  • Identification of edges that appear in all possible MSTs (bridges in weight-buckets)
  • Determination of edges that could appear in at least one MST
  • Computation of each edge's best replacement if it is removed (for alternate MSTs)
  • Calculation of the second-best overall spanning tree cost

Key Components:
  1. Edge and DSU Structures
    • `Edge {u, v, w, id}` holds endpoints, weight, and original index
    • `DSU` supports `find`, `unite`, and specialized linking for path-based replacement queries

  2. Kruskal Sweep with Bucketed Weights
    • Sort all edges by weight and process equal-weight groups (buckets)
    • For each bucket:
      - Build a component-level auxiliary graph from current DSU state
      - Run Tarjan's DFS to detect bridges, these mark edges always in every MST
      - Perform DSU unions to incorporate bucket edges into the MST, recording which edges join
    • Accumulates `mst_weight` and populates the adjacency list `adj` for the resulting MST

  3. LCA Preprocessing for Max-Edge Queries
    • Root MST at node 1 and run DFS to record `depth`, immediate `parent`, and `maxWeight`/`maxEdge` for binary lifting
    • Precomputes `parent[k][v]` and `maxWeight[k][v]` tables in O(n·log(n)) to answer "maximum-edge-on-path" queries in O(log(n))

  4. Edge-by-Edge Replacement Analysis
    • For each original edge not in the MST:
      - Use LCA-based max-edge query on the path between its endpoints to find the heaviest edge in the MST path
      - Compute candidate cost = `mst_weight + w(new) - w(maxOnPath)`
      - If `candidate == mst_weight`, the new edge can swap to yield another MST (`could_be_in_MST`)
      - Otherwise, update the global `second_weight` as the minimum candidate cost > `mst_weight`
    • Applies replacements to record, for each MST edge, the best non-MST edge that could replace it

  5. Best Replacement
    • best_replace[i] = the best edge replacement for edge `i` after removing it from the MST containng `i`
      - equal to -1 if we can not construct any spanning tree without edge `i`

  6. Results and Complexity
    • Arrays:
      - `always_in_MST[id]`, `could_be_in_MST[id]`, `in_mst[id]`, `best_replace[id]` hold classification data per edge
    • Computes both the unique/alternate MST edges and the second-best MST weight
    • Time Complexity: O(m·log(m) + m·α(n) + n·log(n) + m·log(n)) ≈ O(m·log(m))
============================================================================================================*/

struct MST {
    struct Edge {
        int u, v, id; ll w;
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
        void link_child_to_parent(int child, int parent) {
            child = find(child);
            parent = find(parent);
            par[child] = parent;
        }
    };
    
    DSU *dsu_help;
    vector<vector<ll>> maxWeight;
    vector<Edge> edges, edges_org;
    vector<vector<tuple<int, ll, int>>> adj;
    vector<vector<int>> parent, maxEdge;
    vector<int> depth, best_replace, edge_to_parent;
    vector<bool> in_mst, could_be_in_MST, always_in_MST;
    int n, m, LOG; ll mst_weight = 0, second_weight = LLONG_MAX;
    
    void build_mst_and_classify(DSU &dsu) {
        sort(edges.begin(), edges.end(), [] (auto &a, auto &b) {return a.w < b.w;});
        int i = 0;
        while (i < m) {
            int j = i, k = 0; ll w = edges[i].w;
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
        maxWeight.assign(LOG, vector<ll>(n + 1, 0));
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
    void apply_replacement(int u, int anc, int cand_id) {
        while (true) {
            u = dsu_help->find(u);
            if (depth[u] <= depth[anc]) break;
            int id = edge_to_parent[u];
            if (best_replace[id] == -1 or edges_org[best_replace[id]].w > edges_org[cand_id].w) 
                best_replace[id] = cand_id;
            dsu_help->link_child_to_parent(u, parent[0][u]);
        }
    }
    pair<ll, int> get_max_on_path(int u, int v) {
        if (depth[u] < depth[v]) swap(u, v);
        ll best_w = 0; int best_id = -1;
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
            best_replace[e.id] = pr.second;
            ll cand = mst_weight + e.w - pr.first;
            if (cand == mst_weight) could_be_in_MST[e.id] = true;
            else if (cand > mst_weight) second_weight = min(second_weight, cand);
            int u = e.u, v = e.v;
            if (depth[u] < depth[v]) swap(u, v);
            for (int k = LOG - 1; k >= 0; k--) if (depth[u] - (1 << k) >= depth[v]) u = parent[k][u];
            if (u != v) {
                for (int k = LOG - 1; k >= 0; k--) 
                    if (parent[k][u] != parent[k][v]) u = parent[k][u], v = parent[k][v];
                u = parent[0][u];    
            }        
            apply_replacement(e.u, u, e.id);
            apply_replacement(e.v, u, e.id);
        }
    }
    void init() {
        cin >> n >> m;
        edges.resize(m), best_replace.assign(m, -1), edge_to_parent.assign(n + 1, -1);;
        for (int i = 0; i < m; i++) {
            cin >> edges[i].u >> edges[i].v >> edges[i].w;
            edges[i].id = i;
        }
        edges_org = edges;
        in_mst.assign(m, false), could_be_in_MST.assign(m, false), always_in_MST.assign(m, false), adj.assign(n + 1, {});
        DSU dsu(n);
        build_mst_and_classify(dsu);
        preprocess_lca();
        for (int u = 2; u <= n; u++) edge_to_parent[u] = maxEdge[0][u];
        dsu_help = new DSU(n);
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

struct Kirchhoff {
    int n, m;
    vector<vector<int>> A;
    vector<vector<ld>> L, M;

    ll solve() {
        cin >> n >> m;
        A.assign(n, vector<int>(n, 0));
        L.assign(n, vector<ld>(n, 0));
        M.assign(n - 1, vector<ld>(n - 1));
        for (int i = 0, u, v; i < m; i++) {
            cin >> u >> v; u--, v--;
            A[u][v]++, A[v][u]++;
        }
        for (int i = 0; i < n; i++) {
            int deg = 0;
            for (int j = 0; j < n; j++) deg += A[i][j];
            L[i][i] = (ld) deg;
            for (int j = 0; j < n; j++) if (i != j)
                L[i][j] = -(ld) A[i][j];
        }
        for (int i = 0; i < n - 1; i++)
            for (int j = 0; j < n - 1; j++) M[i][j] = L[i][j];
        ld det = 1;
        for (int i = 0; i < n - 1; i++) {
            int pivot = i;
            for (int r = i + 1; r < n - 1; r++) 
                if (fabsl(M[r][i]) > fabsl(M[pivot][i])) pivot = r;
            if (fabsl(M[pivot][i]) < 1e-15L) {det = 0; break;}
            if (pivot != i) {swap(M[pivot], M[i]); det = -det;}
            det *= M[i][i];
            ld inv = 1.0L / M[i][i];
            for (int r = i + 1; r < n - 1; r++) {
                ld factor = M[r][i] * inv;
                for (int c = i; c < n - 1; c++) 
                    M[r][c] -= factor * M[i][c];
            }
        } 
        ll ans = (ll) (det + (det > 0 ? 0.5L : -0.5L));
        return ans;
    }
};
