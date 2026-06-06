#include <bits/stdc++.h>
using namespace std;
using ll = long long int;

/*============================================================================================================
Graph Connectivity & Decomposition Toolkit (Undirected + Directed)

What this file contains
  • A collection of classic graph algorithms related to connectivity, critical elements,
    and decomposition, written in competitive-programming style.
  • Covers both static and dynamic settings, undirected and directed graphs.

Included Algorithms:

1) Connected Components (Undirected):
  • Finds all maximal connected subgraphs using iterative DFS
  • Useful as a basic preprocessing step
  • Time: O(n + m)

2) Bridges (Undirected, Offline):
  • Detects edges whose removal increases the number of connected components
  • Uses DFS discovery time + low-link values
  • Handles multi-edges correctly by skipping the parent edge once
  • Time: O(n + m)

3) Bridges (Undirected, Online / Incremental):
  • Maintains the number of bridges while edges are added dynamically
  • Uses DSU over 2-edge-connected components + parent pointers + LCA-style merging
  • Much harder than offline version, but supports online updates
  • Time: O((n + m)·log(n))

4) Articulation Points (Cut-Vertices):
  • Vertices whose removal increases the number of connected components
  • DFS-based, closely related to bridge detection
  • Root rule: root is articulation if it has >1 DFS children
  • Non-root rule: u is articulation if ∃ child v with low[v] ≥ in[u]
  • Time: O(n + m)

5) Strongly Connected Components (Directed, Kosaraju):
  • Decomposes a directed graph into SCCs
  • Two DFS passes: order on original graph, then DFS on reversed graph
  • Produces SCC list and condensation DAG
  • Time: O(n + m)

6) Strong Orientation (Undirected → Directed):
  • Assigns a direction to each undirected edge
  • Goal: minimize number of SCCs after orientation
  • Key fact: number of SCCs = connected components + number of bridges
  • Bridges force direction, non-bridge edges can be oriented arbitrarily
  • Time: O(n + m)

7) Edge Connectivity (λ) & Vertex Connectivity (κ):
  • Edge connectivity λ:
      minimum number of edges whose removal disconnects the graph
  • Vertex connectivity κ:
      minimum number of vertices whose removal disconnects the graph
  • Key inequalities (Whitney):
      κ ≤ λ ≤ δ   (δ = minimum degree)
  • Approaches:
    a) Pairwise max-flow (conceptual / brute-force)
       – Compute min s-t cut for all pairs
       – Node-splitting trick for vertex connectivity
    b) Stoer–Wagner global min-cut
       – Finds λ directly without enumerating pairs
       – Works on weighted undirected graphs
       – Time: O(V³) or O(V·E)

General Notes
  • All graphs are 0-indexed unless stated otherwise
  • All DFS-based algorithms are linear in n + m
  • Online bridge algorithm is significantly more complex but supports incremental updates
  • SCC condensation graph is a DAG and may contain duplicate edges unless filtered

Typical Use-Cases
  • Network reliability and fault detection
  • Dynamic graph monitoring
  • Graph compression and decomposition
  • Preprocessing for higher-level algorithms (flows, orientations, connectivity proofs)
============================================================================================================*/

struct FindCC {
    int n;
    vector<int> comp;
    vector<bool> used;
    vector<vector<int>> G;

    void DFS(const int u) {
        stack<int> st;
        st.push(u);
        while (!st.empty()) {
            int cur = st.top();
            st.pop();
            if (!used[cur]) {
                used[cur] = true;
                comp.push_back(cur);
                for (int i = static_cast<int>(G[cur].size()) - 1; i >= 0; i--) st.push(G[cur][i]);
            }
        }
    }
    void find_cc() {
        used.assign(n, false);
        for (int u = 0; u < n; u++) if (!used[u]) { comp.clear(); DFS(u); }
    }
};

struct FindBridges {
    int n, timer;
    vector<bool> vis;
    vector<int> in, low;
    vector<vector<int>> G;

    void is_bridge(int u, int v);

    void DFS(const int u, const int p = -1) {
        vis[u] = true, in[u] = low[u] = timer++;
        bool parent_skipped = false;
        for (const int v : G[u]) {
            if (v == p and !parent_skipped) { parent_skipped = true; continue; }
            if (vis[v]) low[u] = min(low[u], in[v]);
            else {
                DFS(v, u);
                low[u] = min(low[u], low[v]);
                if (low[v] > in[u]) is_bridge(u, v);
            }
        }
    }
    void find_bridges() {
        timer = 0, vis.assign(n, false), in.assign(n, -1), low.assign(n, -1);
        for (int u = 0; u < n; u++) if (!vis[u]) DFS(u);
    }
};

struct FindBridgesOnline {
    int bridges, lca_iteration;
    vector<int> par, dsu_2ecc, dsu_cc, dsu_cc_size, last_visit;

    void init(const int n) {
        par.resize(n), dsu_2ecc.resize(n), dsu_cc.resize(n);
        dsu_cc_size.resize(n), last_visit.assign(n, 0), lca_iteration = bridges = 0;
        for (int i = 0; i < n; i++)
            dsu_2ecc[i] = i, dsu_cc[i] = i, dsu_cc_size[i] = 1, par[i] = -1;
    }
    int find_2ecc(const int u) {
        if (u == -1) return -1;
        return dsu_2ecc[u] == u ? u : dsu_2ecc[u] = find_2ecc(dsu_2ecc[u]);
    }
    int find_cc(int u) {
        u = find_2ecc(u);
        return dsu_cc[u] == u ? u : dsu_cc[u] = find_cc(dsu_cc[u]);
    }
    void make_root(int u) {
        const int root = u;
        int child = -1;
        while (u != -1) {
            const int p = find_2ecc(par[u]);
            par[u] = child, dsu_cc[u] = root, child = u, u = p;
        }
        dsu_cc_size[root] = dsu_cc_size[child];
    }
    void merge_path(int a, int b) {
        lca_iteration++;
        vector<int> path_a, path_b;
        int lca = -1;
        while (lca == -1) {
            if (a != -1) {
                a = find_2ecc(a), path_a.push_back(a);
                if (last_visit[a] == lca_iteration) { lca = a; break; }
                last_visit[a] = lca_iteration, a = par[a];
            }
            if (b != -1) {
                b = find_2ecc(b), path_b.push_back(b);
                if (last_visit[b] == lca_iteration) { lca = b; break; }
                last_visit[b] = lca_iteration, b = par[b];
            }
        }
        for (const int u : path_a) {
            dsu_2ecc[u] = lca;
            if (u == lca) break;
            bridges--;
        }
        for (const int u : path_b) {
            dsu_2ecc[u] = lca;
            if (u == lca) break;
            bridges--;
        }
    }
    void add_edge(int a, int b) {
        a = find_2ecc(a), b = find_2ecc(b);
        if (a == b) return;
        int ca = find_cc(a);
        if (int cb = find_cc(b); ca != cb) {
            bridges++;
            if (dsu_cc_size[ca] > dsu_cc_size[cb]) swap(a, b), swap(ca, cb);
            make_root(a);
            par[a] = dsu_cc[a] = b;
            dsu_cc_size[cb] += dsu_cc_size[a];
        }
        else merge_path(a, b);
    }
};

struct FindCutpoints {
    int n, timer;
    vector<bool> vis;
    vector<int> in, low;
    vector<vector<int>> G;

    void is_cutpoint(int u);

    void DFS(const int u, const int p = -1) {
        int child = 0;
        vis[u] = true, in[u] = low[u] = timer++;
        for (const int v : G[u]) {
            if (v == p) continue;
            if (vis[v]) low[u] = min(low[u], in[v]);
            else {
                DFS(v, u);
                low[u] = min(low[u], low[v]);
                if (low[v] >= in[u] and p != -1) is_cutpoint(u);
                child++;
            }
        }
        if (p == -1 and child > 1) is_cutpoint(u);
    }
    void find_cutpoints() {
        timer = 0, vis.assign(n, false), in.assign(n, -1), low.assign(n, -1);
        for (int u = 0; u < n; u++) if (!vis[u]) DFS(u);
    }
};

struct FindSCC {
    int n;
    vector<bool> vis;
    vector<vector<int>> G, G_cond, comps;

    void DFS(const int u, const vector<vector<int>> &G, vector<int> &order){
        vis[u] = true;
        for (const int v : G[u]) if (!vis[v]) DFS(v, G, order);
        order.push_back(u);
    }
    void find_scc() {
        vector<int> order;
        vis.assign(n, false);
        for (int u = 0; u < n; u++) if (!vis[u]) DFS(u, G, order);
        vector<vector<int>> G_rev(n);
        for (int u = 0; u < n; u++) for (const int v : G[u]) G_rev[v].push_back(u);
        vis.assign(n, false);
        reverse(order.begin(), order.end());
        vector roots(n, 0);
        for (const int u : order)
            if (!vis[u]) {
                vector<int> comp;
                DFS(u, G_rev, comp);
                comps.push_back(comp);
                const int root = *min_element(comp.begin(), comp.end());
                for (const int v : comp) roots[v] = root;
            }
        G_cond.assign(n, {});
        for (int u = 0; u < n; u++) for (const int v : G[u]) if (roots[u] != roots[v]) G_cond[roots[u]].push_back(roots[v]);
    }
};

struct StrongOrientation {
    int n, m, bridges;
    vector<char> orientation;
    vector<int> in, low;
    vector<bool> edge_used;
    vector<pair<int, int>> edges;
    vector<vector<pair<int, int>>> G;

    void find_bridges(const int u) {
        static int timer = 0;
        low[u] = in[u] = timer++;
        for (const auto [fst, snd] : G[u]) {
            if (edge_used[snd]) continue;
            edge_used[snd] = true;
            orientation[snd] = u == edges[snd].first ? '>' : '<';
            if (in[fst] == -1) {
                find_bridges(fst);
                low[u] = min(low[u], low[fst]);
                if (low[fst] > in[u]) bridges++;
            }
            else low[u] = min(low[u], in[fst]);
        }
    }
    void init() {
        bridges = 0, G.assign(n, {}), edges.resize(m), in.assign(n, -1);
        low.assign(n, -1), orientation.resize(m), edge_used.assign(m, false);
        for (int i = 0, u, v; i < m; i++) {
            cin >> u >> v;
            G[--u].emplace_back(--v, i);
            G[v].emplace_back(u, i);
            edges[i] = {u, v};
        }
        int comps = 0;
        for (int u = 0; u < n; u++) if (in[u] == -1) { comps++; find_bridges(u); }
    }
};

struct StoerWagner {
    int n, m;
    vector<int> vert;
    vector<vector<ll>> weight;

    ll solve() {
        vert.resize(n);
        ll best = LLONG_MAX;
        iota(vert.begin(), vert.end(), 0);
        for (int phase = n; phase > 1; phase--) {
            int prev = 0;
            vector<ll> w(phase, 0);
            vector added(phase, false);
            for (int i = 0; i < phase; i++) {
                int sel = -1;
                for (int j = 0; j < phase; j++)
                    if (!added[j] and (sel == -1 or w[j] > w[sel])) sel = j;
                added[sel] = true;
                if (i == phase - 1) {
                    best = min(best, w[sel]);
                    for (int j = 0; j < phase; j++) {
                        weight[vert[prev]][vert[j]] += weight[vert[sel]][vert[j]];
                        weight[vert[j]][vert[prev]] = weight[vert[prev]][vert[j]];
                    }
                    vert.erase(vert.begin() + sel);
                } else {
                    prev = sel;
                    for (int j = 0; j < phase; j++)
                        if (!added[j]) w[j] += weight[vert[sel]][vert[j]];
                }
            }
        }
        return best;
    }
    void input() {
        cin >> n >> m;
        weight.assign(n, vector<ll>(n, 0));
        for (int i = 0; i < m; i++) {
            int u, v; ll w;
            cin >> u >> v >> w; u--, v--;
            weight[u][v] += w, weight[v][u] += w;
        }
        solve();
    }
};
