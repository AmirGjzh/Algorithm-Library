#include <bits/stdc++.h>
using namespace std;
using ll = long long int;

/*============================================================================================================
Connected Components (undirected)

Description:
  • Partition an undirected graph into its maximal connected subgraphs
  • Here we use an explicit stack (iterative DFS)

Applications:
  • Determining cluster membership
  • Preprocessing for bridge or articulation‑point algorithms
  • Counting connected pieces in grid or mesh based problems. 

Notes:
  • Each time we find an unvisited u, we clear the current component vector and DFS to collect all its vertices
  • We push neighbors in reverse order so that smaller-index neighbors are processed first (if you care about lexicographic order)

Order: O(n + m)
============================================================================================================*/

struct FindCC {
    int n;
    vector<int> comp;
    vector<bool> used;
    vector<vector<int>> G;

    void DFS(int u) {
        stack<int> st;
        st.push(u);
        while (st.size()) {
            int cur = st.top();
            st.pop();
            if (!used[cur]) {
                used[cur] = true;
                comp.push_back(cur);
                for (int i = G[cur].size() - 1; i >= 0; i--) st.push(G[cur][i]);
            }
        }
    }
    void find_cc() {
        used.assign(n, false);
        for (int u = 0; u < n; u++) if (!used[u]) {comp.clear(); DFS(u);}
    }
};

/*============================================================================================================
Bridges in an undirected graph

Description:
  • An edge is a bridge if its removal increases the number of connected components
  • Uses DFS with discovery times `in[u]` and low‑link values `low[u]`

Applications:
  • Network reliability: identify single points of failure
  • Graph augmentation: guaranteed edges to add for 2‑edge‑connectivity
  • Preprocessing for strongly oriented or biconnected‑component algorithms

Notes:
  • Must ignore the parent edge when scanning adjacency to avoid treating it as a back‑edge
  • In the presence of parallel edges, you must skip only one occurrence of the parent—extra copies become non‑bridges

Order: O(n + m)
============================================================================================================*/

struct FindBridges {
    int n, timer;
    vector<bool> vis;
    vector<int> in, low;
    vector<vector<int>> G;

    void is_bridge(int u, int v);

    void DFS(int u, int p = -1) {
        vis[u] = true, in[u] = low[u] = timer++;
        bool parent_skipped = false;
        for (int v : G[u]) {
            if (v == p and !parent_skipped) {parent_skipped = true; continue;}
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

/*============================================================================================================
Bridges (Online)

Description:
  • Maintains the bridge count as edges are added one by one
  • Uses DSU over “2‑eccentric components” + parent pointers + LCA‑style merging

Applications:
  • Dynamic connectivity in incremental networks
  • Real‑time monitoring of network vulnerability

Notes:
  • `init(n)` resets DSU arrays and counters
  • `add_edge(a, b)` finds the 2‑ecc representatives, merges CCs or merges paths via LCA logic

Order: O((n + m)·log(n))
============================================================================================================*/

struct FindBridgesOnline {
    int bridges, lca_iteration;
    vector<int> par, dsu_2ecc, dsu_cc, dsu_cc_size, last_visit;

    void init(int n) {
        par.resize(n), dsu_2ecc.resize(n), dsu_cc.resize(n);
        dsu_cc_size.resize(n), last_visit.assign(n, 0), lca_iteration = bridges = 0;
        for (int i = 0; i < n; i++) 
            dsu_2ecc[i] = i, dsu_cc[i] = i, dsu_cc_size[i] = 1, par[i] = -1;
    }
    int find_2ecc(int u) {
        if (u == -1) return -1;
        return dsu_2ecc[u] == u ? u : dsu_2ecc[u] = find_2ecc(dsu_2ecc[u]);    
    }
    int find_cc(int u) {
        u = find_2ecc(u);
        return dsu_cc[u] == u ? u : dsu_cc[u] = find_cc(dsu_cc[u]);
    }
    void make_root(int u) {
        int root = u, child = -1;
        while (u != -1) {
            int p = find_2ecc(par[u]);
            par[u] = child, dsu_cc[u] = root, child = u, u = p;
        }
        dsu_cc_size[root] = dsu_cc_size[child];
    }
    void merge_path(int a, int b) {
        lca_iteration++;
        vector<int> path_a, path_b;
        int lca = -1;
        while (lca == -1 ) {
            if (a != -1) {
                a = find_2ecc(a), path_a.push_back(a);
                if (last_visit[a] == lca_iteration) {lca = a; break;}
                last_visit[a] = lca_iteration, a = par[a];
            }
            if (b != -1) {
                b = find_2ecc(b), path_b.push_back(b);
                if (last_visit[b] == lca_iteration) {lca = b; break;}
                last_visit[b] = lca_iteration, b = par[b];
            }
        }
        for (int u : path_a) {
            dsu_2ecc[u] = lca;
            if (u == lca) break;
            bridges--;
        }
        for (int u : path_b) {
            dsu_2ecc[u] = lca;
            if (u == lca) break;
            bridges--;
        }
    }
    void add_edge(int a, int b) {
        a = find_2ecc(a), b = find_2ecc(b);
        if (a == b) return;
        int ca = find_cc(a), cb = find_cc(b);
        if (ca != cb) {
            bridges++;
            if (dsu_cc_size[ca] > dsu_cc_size[cb]) swap(a, b), swap(ca, cb);
            make_root(a);
            par[a] = dsu_cc[a] = b;
            dsu_cc_size[cb] += dsu_cc_size[a];
        }
        else merge_path(a, b);
    }
};

/*============================================================================================================
Cut‑vertices (Articulation Points)

Description:
  • A vertex whose removal increases the number of connected components
  • Very similar to bridge‑finding: track child count + low‑link vs entry time

Applications:
  • Identifying critical routers or junctions
  • Graph biconnectivity decomposition
  • Preprocessing for 2‑vertex‑connected augmentation

Notes:
  • For root u (p == –1), u is a cut‑point if it has more than one DFS child
  • For non‑root u, u is a cut‑point whenever there exists a child v with low[v] ≥ in[u]
  • Important: for some u, `is_cutpoint(u)` may be invoked multiple times—once per qualifying child v. 

Order: O(n + m)
============================================================================================================*/

struct FindCutpoints {
    int n, timer;
    vector<bool> vis;
    vector<int> in, low;
    vector<vector<int>> G;

    void is_cutpoint(int u);

    void DFS(int u, int p = -1) {
        int child = 0;
        vis[u] = true, in[u] = low[u] = timer++;
        for (int v : G[u]) {
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

/*============================================================================================================
Strongly Connected Components (Kosaraju)

Description:
  • Decompose a directed graph into maximal strongly connected subgraphs (SCCs)
  • Two‑pass DFS: first pass to compute exit order, second pass on reversed graph to extract components

Applications:
  1. Deadlock detection
  2. Component‑level DAG scheduling
  3. Graph compression / meta‑graph construction

Notes:
  • Use separate visited arrays for each pass to avoid confusion
  • After collecting `comp`, choose a unique representative `root = min(comp)` before building the condensation graph
  • When adding edges to `g_cond`, check for and skip duplicates to keep the DAG clean

Order: O(n + m)
============================================================================================================*/

struct FindSCC {
    int n;
    vector<bool> vis;
    vector<vector<int>> G, G_cond, comps;

    void DFS(int u, const vector<vector<int>> &G, vector<int> &order){
        vis[u] = true;
        for (int v : G[u]) if (!vis[v]) DFS(v, G, order);
        order.push_back(u);        
    }
    void find_scc() {
        vector<int> order;
        vis.assign(n, false);
        for (int u = 0; u < n; u++) if (!vis[u]) DFS(u, G, order);
        vector<vector<int>> G_rev(n);
        for (int u = 0; u < n; u++) for (int v : G[u]) G_rev[v].push_back(u);
        vis.assign(n, false);
        reverse(order.begin(), order.end());
        vector<int> roots(n, 0);
        for (int u : order) 
            if (!vis[u]) {
                vector<int> comp;
                DFS(u, G_rev, comp);
                comps.push_back(comp);
                int root = *min_element(comp.begin(), comp.end());
                for (int v : comp) roots[v] = root; 
            }
        G_cond.assign(n, {});
        for (int u = 0; u < n; u++) for (int v : G[u]) if (roots[u] != roots[v]) G_cond[roots[u]].push_back(roots[v]);
    }
};

/*============================================================================================================
Strong Orientation

Description:
  • Assigns a direction (‘>’ or ‘<’) to each undirected edge so the resulting digraph
    has as few SCCs as possible (ideal: 1 SCC if no bridges)
  • First identify bridges while orienting edges along a DFS tree, non‑bridge edges can be oriented arbitrarily

Applications:
  1. Making communication networks strongly connected
  2. Generating a minimal‑feedback orientation
  3. Designing one‑way street systems with minimal disconnections

Notes:
  • The number of SCCs after orientation = CCs + bridges
  • ‘>’ means the direction is the same as input, ‘<’ means the direction is reversed
  • `init()` reads edges, builds adjacency with edge IDs, and resets all arrays
  • `dfs_orient_and_count_bridges(u)` both orients edges and increments `bridges` when `low[v] > in[u]`
  • After DFS, any edge not yet in `orientation` can be assigned arbitrarily

Order: O(n + m)
============================================================================================================*/

struct StrongOrientation {
    int n, m, bridges;
    vector<char> orientation;
    vector<int> in, low;
    vector<bool> edge_used;
    vector<pair<int, int>> edges;
    vector<vector<pair<int, int>>> G;

    void find_bridges(int u) {
        static int timer = 0;
        low[u] = in[u] = timer++;
        for (auto p : G[u]) {
            if (edge_used[p.second]) continue;
            edge_used[p.second] = true;
            orientation[p.second] = (u == edges[p.second].first ? '>' : '<');
            if (in[p.first] == -1) {
                find_bridges(p.first);
                low[u] = min(low[u], low[p.first]);
                if (low[p.first] > in[u]) bridges++;
            }
            else low[u] = min(low[u], in[p.first]);
        }
    }
    void init() {
        bridges = 0, G.assign(n, {}), edges.resize(m), in.assign(n, -1);
        low.assign(n, -1), orientation.resize(m), edge_used.assign(m, false);
        for (int i = 0, u, v; i < m; i++) {
            cin >> u >> v;
            G[--u].push_back({--v, i});
            G[v].push_back({u, i});
            edges[i] = {u, v};
        }
        int comps = 0;
        for (int u = 0; u < n; u++) if (in[u] == -1) {comps++; find_bridges(u);}
    }
};

/*============================================================================================================
Edge & Vertex Connectivity on undirected graphs

Definitions:
  • Edge connectivity λ: the smallest number of edges whose removal makes G disconnected
    – Equivalently, λ = minₛₜ (size of a minimum s–t edge cut)
  • Vertex connectivity κ: the smallest number of vertices whose removal makes G disconnected
    – Equivalently, κ = minₛₜ (size of a minimum s–t vertex cut)

Key facts:
  • Whitney’s inequalities:  
    κ ≤ λ ≤ δ  (δ = min vertex degree)
  • Ford–Fulkerson tells us:  
    maximum number of edge‑disjoint s–t paths = min s–t edge cut

Implementation sketches:
  1. Brute-force via max‑flow
    - For every ordered pair (s,t):  
      • Build a flow network with unit capacities on each original edge
      • Compute max‐flow(s→t) (e.g. Edmonds‑Karp or preferably Dinic) 
      • Track the minimum over all pairs → this is λ
    - To get κ, first **node‑split** every non {s,t} vertex x into x_in → x_out with capacity=1,  
      redirect edges through these, then run the same pairwise max‐flow

  2. Global min‐cut (Stoer–Wagner)  
    - Finds λ in O(V³) or O(VE) without checking all pairs
    - Repeatedly “merge” the most tightly connected vertices, tracking the cut weight

Practical tips:
  • Use Dinic (O(E.√V)) or Push–Relabel for faster flows on dense graphs
  • For vertex cuts, remember to:
    – Create two nodes per original node
    – Give the splitting edge capacity = 1
    – Redirect each original (u,v) into u_out → v_in and v_out → u_in  
  • Stoer–Wagner is simple to code and often beats pairwise flows on large sparse graphs

Complexities:
  • Pairwise max‑flow: O(V²·FlowTime) (worst‑case O(V³E²) with Edmonds‑Karp)
  • Stoer–Wagner: O(V³) or O(V.E)  
  • Node‑splitting doubles V, so flows run on ≈2V vertices
============================================================================================================*/

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
            vector<bool> added(phase, false);
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