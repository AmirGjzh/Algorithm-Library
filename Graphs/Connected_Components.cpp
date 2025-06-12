#include <bits/stdc++.h>
using namespace std;

/*--------------------------------------------------------------------------------------------------------
Connected Components in undirected graph:
simply using DFS or BFS
this code implent DFS with stack
Order = O(n + m)
--------------------------------------------------------------------------------------------------------*/

struct FindCC {
    int n;
    vector<int> comp;
    vector<bool> used;
    vector<vector<int>> g;

    void DFS(int u) {
        stack<int> st;
        st.push(u);
        while (st.size()) {
            int cur = st.top();
            st.pop();
            if (!used[cur]) {
                used[cur] = true;
                comp.push_back(cur);
                for (int i = g[cur].size() - 1; i >= 0; i--)
                    st.push(g[cur][i]);
            }
        }
    }

    void find_cc() {
        used.assign(n, false);
        for (int u = 0; u < n; u++) 
            if (!used[u]) {
                comp.clear();
                DFS(u);
            }
    }
};

/*--------------------------------------------------------------------------------------------------------
Bridges:
Bridges are the edges that whan we remove it, number of CCs will increase
We run DFS to find all bridges:
edge (u, to) is a bridge if and only if any backedges from to (or its descendants) goes to an ancestor of u
note that in case of multiple edges, we must ignore the parent edge once!
Order = O(n + m)
--------------------------------------------------------------------------------------------------------*/

struct FindBridges {
    int n, timer;
    vector<bool> vis;
    vector<int> in, low;
    vector<vector<int>> g;

    void is_bridge(int u, int v);

    void DFS(int u, int p = -1) {
        vis[u] = true;
        in[u] = low[u] = timer++;
        bool parent_skipped = false;
        for (int v : g[u]) {
            if (v == p and !parent_skipped) {
                parent_skipped = true;
                continue;
            }
            if (vis[v]) low[u] = min(low[u], in[v]);
            else {
                DFS(v, u);
                low[u] = min(low[u], low[v]);
                if (low[v] > in[u]) is_bridge(u, v);
            }
        }
    }

    void find_bridges() {
        timer = 0;
        vis.assign(n, false);
        in.assign(n, -1);
        low.assign(n, -1);
        for (int u = 0; u < n; u++) 
            if (!vis[u])
                DFS(u);
    }
};

/*--------------------------------------------------------------------------------------------------------
Finding Bridges Online:
This was hard, 
So just the implementation and how it works :)
just call init and then add edges 
Order = O((n + m).Log(n))
--------------------------------------------------------------------------------------------------------*/

struct FindBridgesOnline {
    int bridges, lca_iteration;
    vector<int> par, dsu_2ecc, dsu_cc, dsu_cc_size, last_visit;

    void init(int n) {
        par.resize(n);
        dsu_2ecc.resize(n);
        dsu_cc.resize(n);
        dsu_cc_size.resize(n);
        last_visit.assign(n, 0);
        lca_iteration = bridges = 0;
        for (int i = 0; i < n; i++) {
            dsu_2ecc[i] = i;
            dsu_cc[i] = i;
            dsu_cc_size[i] = 1;
            par[i] = -1;
        }
    }

    int find_2ecc(int u) {
        if (u == -1)
            return -1;
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
            par[u] = child;
            dsu_cc[u] = root;
            child = u;
            u = p;
        }
        dsu_cc_size[root] = dsu_cc_size[child];
    }

    void merge_path(int a, int b) {
        lca_iteration++;
        vector<int> path_a, path_b;
        int lca = -1;
        while (lca == -1 ) {
            if (a != -1) {
                a = find_2ecc(a);
                path_a.push_back(a);
                if (last_visit[a] == lca_iteration) {
                    lca = a;
                    break;
                }
                last_visit[a] = lca_iteration;
                a = par[a];
            }
            if (b != -1) {
                b = find_2ecc(b);
                path_b.push_back(b);
                if (last_visit[b] == lca_iteration) {
                    lca = b;
                    break;
                }
                last_visit[b] = lca_iteration;
                b = par[b];
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
        a = find_2ecc(a);
        b = find_2ecc(b);
        if (a == b) return;
        int ca = find_cc(a);
        int cb = find_cc(b);
        if (ca != cb) {
            bridges++;
            if (dsu_cc_size[ca] > dsu_cc_size[cb]) {
                swap(a, b);
                swap(ca, cb);
            }
            make_root(a);
            par[a] = dsu_cc[a] = b;
            dsu_cc_size[cb] += dsu_cc_size[a];
        }
        else merge_path(a, b);
    }
};

/*--------------------------------------------------------------------------------------------------------
Cut points or vertices:
these vertices are the on that if we remove them, number of CCs will increase, the code is veryyy similar
to FindBridges
Note: for sum u the code may call is_cutpoint(u) many times (for some of u's children)
Order = O(n + m)
--------------------------------------------------------------------------------------------------------*/

struct FindCutpoints {
    int n, timer;
    vector<bool> vis;
    vector<int> in, low;
    vector<vector<int>> g;

    void is_cutpoint(int u);

    void DFS(int u, int p = -1) {
        int child = 0;
        vis[u] = true;
        in[u] = low[u] = timer++;
        for (int v : g[u]) {
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
        timer = 0;
        vis.assign(n, false);
        in.assign(n, -1);
        low.assign(n, -1);
        for (int u = 0; u < n; u++) 
            if (!vis[u])
                DFS(u);
    }
};

/*--------------------------------------------------------------------------------------------------------
Strongly Connected Components:
Finding SCC in a directed graph 
first we run a dfs to obtain the edges in desc order of their exit time from dfs
then run dfs on the=at order on the transpose graph
we also create the condensation graph witch each vertex in it is a SCC and it is a DAG
Order = O(n + m)
--------------------------------------------------------------------------------------------------------*/

struct FindSCC {
    int n;
    vector<bool> vis;
    vector<vector<int>> g, g_cond, comps;

    void DFS(int u, vector<vector<int>> &g, vector<int> &order){
        vis[u] = true;
        for (int v : g[u])
            if (!vis[v])
                DFS(v, g, order);
        order.push_back(u);        
    }

    void find_scc() {
        vector<int> order;
        vis.assign(n, false);
        for (int u = 0; u < n; u++)
            if (!vis[u])
                DFS(u, g, order);
        vector<vector<int>> g_rev(n);
        for (int u = 0; u < n; u++)
            for (int v : g[u])
                g_rev[v].push_back(u);
        vis.assign(n, false);
        reverse(order.begin(), order.end());
        vector<int> roots(n, 0);
        for (int u : order) {
            if (!vis[u]) {
                vector<int> comp;
                DFS(u, g_rev, comp);
                comps.push_back(comp);
                int root = *min_element(comp.begin(), comp.end());
                for (int v : comp) roots[v] = u; 
            }
        }
        g_cond.assign(n, {});
        for (int u = 0; u < n; u++)
            for (int v : g[u])
                if (roots[u] != roots[v])
                    g_cond[roots[u]].push_back(roots[v]);
    }
};

/*--------------------------------------------------------------------------------------------------------
Strongly Orientation:
assigning a direction to each edge of a graph to make it strongly connected
but, we can do it if the graph does not have any bridges
in general we can make the number of SCCs minimal by removing bridges and what left is some CCs that can be
oriented
the number is = CCs + bridges
> means the direction is similar to input and < meanse reverse
Order = O(n + m)?
--------------------------------------------------------------------------------------------------------*/

struct StrongOrientation {
    int n, m, bridges;
    string orientation;
    vector<int> in, low;
    vector<bool> edge_used;
    vector<pair<int, int>> edges;
    vector<vector<pair<int, int>>> g;

    void find_bridges(int u) {
        static int timer = 0;
        low[u] = in[u] = timer++;
        for (auto p : g[u]) {
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
        bridges = 0;
        g.assign(n, {});
        edges.resize(m);
        in.assign(n, -1);
        low.assign(n, -1);
        orientation.resize(m);
        edge_used.assign(m, false);
        for (int i = 0, u, v; i < m; i++) {
            cin >> u >> v;
            g[--u].push_back({--v, i});
            g[v].push_back({u, i});
            edges[i] = {u, v};
        }
        int comps = 0;
        for (int u = 0; u < n; u++)
            if (in[u] == -1) {
                comps++;
                find_bridges(u);
            }
    }
};