#include <bits/stdc++.h>
using namespace std;

/*============================================================================================================
Lowest Common Ancestor (LCA) – Three Common Approaches

Description:
  The Lowest Common Ancestor (LCA) problem finds, for any two nodes u and v in a rooted tree,
  their deepest shared ancestor, Below are three well‑known preprocessing/query strategies:

1) Euler Tour + Segment Tree
  • Build:
    – Perform an Euler tour (DFS) to record a sequence of 2n–1 node visits and depths
    – Record first[u] = index of u’s first appearance in the tour
    – Build a segment tree over the tour array so that each node stores the index of the
      shallower of its two children (min‑by‑depth)
    – Preprocessing time: O(n) for DFS + O(n) to build segment tree (4n space)
  • Query:
    – Query the segment tree for the minimum‑depth node in [l..r] in O(log(n))

2) Euler Tour + Sparse Table (RMQ)
  • Build:
    – Same Euler tour and first[u] recording as above
    – Precompute floor(log₂i) for i up to 2n
    – Build a Sparse Table 
    – Preprocessing time: O(n) for DFS + O(n.log(n)) for Sparse Table
  • Query:
    – For LCA(u,v), set [l,r] = sorted(first[u], first[v])
    – Let k = ⌊log₂(r−l+1)⌋; answer = min‑by‑depth( st[k][l], st[k][r−2^k+1] ) in O(1)

3) Binary Lifting (Ancestor Table + Timestamps)
  • Build:
    – Compute LOG = ⌈log₂N⌉.
    – Run a DFS to assign tin[u]/tout[u] timestamps and record up[u][0] = parent
    – Build up[u][k] = 2^k‑th ancestor by up[u][k−1] lookups: O(n.log(n))
  • Query:
    – If u is ancestor of v (checked via tin/tout), return u (and vice versa)
    – Otherwise, for k from LOG down to 0, lift u = up[u][k] whenever up[u][k] is not ancestor of v
    – Finally return up[u][0]. Total per-query: O(log(n))

Complexity Comparison:
  • Space:      Segment Tree O(n),     Sparse Table O(n.log(n)),  Binary Lifting O(n.log(n))
  • Preprocess: Segment Tree O(n),     Sparse Table O(n.log(n)),  Binary Lifting O(n.log(n))
  • Query:      Segment Tree O(log(n)), Sparse Table O(1),        Binary Lifting O(log(n))

When to Use:
  – Choose Sparse Table for fastest queries when memory/time budget allows O(n.log(n)) build
  – Use Segment Tree if memory is tight (linear space) but O(log(n)) query is fine
  – Use Binary Lifting if you also need ancestor jumps or dynamic depth/parent queries alongside LCA
============================================================================================================*/

struct LCA {
    int n;
    vector<vector<int>> st;
    vector<int> h, euler, first, lg;

    LCA(const vector<vector<int>> &g, int root) {
        n = g.size(), h.resize(n);
        first.resize(n), euler.reserve(n << 1);
        DFS(g, root);
        int m = euler.size(); lg.assign(m + 1, 0);
        for (int i = 2; i <= m; i++) lg[i] = lg[i / 2] + 1;
        st.assign(lg[m] + 1, vector<int>(m));
        for (int i = 0; i < m; i++) st[0][i] = euler[i];
        for (int k = 1; k < lg[m] + 1; k++) {
            int len = (1 << k);
            for (int i = 0; i + len <= m; i++) {
                int x = st[k - 1][i], y = st[k - 1][i + (len >> 1)];
                st[k][i] = (h[x] < h[y] ? x : y);
            }
        }
    }
    void DFS(const vector<vector<int>> &g, int u, int p = -1, int height = 0) {
        h[u] = height, first[u] = euler.size();
        euler.push_back(u);
        for (int v : g[u]) 
            if (v != p) {
                DFS(g, v, u, height + 1);
                euler.push_back(u);
            }
    }
    int query(int l, int r) {
        int k = lg[r - l + 1], x = st[k][l], y = st[k][r - (1 << k) + 1];
        return (h[x] < h[y] ? x : y);
    }
    int lca(int u, int v) {
        int l = first[u], r = first[v];
        if (l > r) swap(l, r);
        return query(l, r);
    }
};

struct LCA {
    int n;
    vector<bool> vis;
    vector<int> h, euler, first, seg;

    LCA(vector<vector<int>> &g, int root) {
        n = g.size(), h.resize(n), first.resize(n);
        euler.reserve(n << 1), vis.assign(n, false);
        DFS(g, root);
        int m = euler.size(); seg.resize(m << 2);
        build(0, m - 1);
    }
    void DFS(vector<vector<int>> &g, int u, int height = 0) {
        vis[u] = true, h[u] = height, first[u] = euler.size();
        euler.push_back(u);
        for (auto v : g[u]) 
            if (!vis[v]) {
                DFS(g, v, height + 1);
                euler.push_back(u);
            }
    }
    void build(int l, int r, int id = 1) {
        if (l == r) seg[id] = euler[l];
        else {
            int mid = (l + r) >> 1;
            build(l, mid, id << 1), build(mid + 1, r, id << 1 | 1);
            int left = seg[id << 1], right = seg[id << 1 | 1];
            seg[id] = (h[left] < h[right] ? left : right);
        }
    }
    int query(int L, int R, int l, int r, int id = 1) {
        if (l > R or r < L) return -1;
        if (l >= L and r <= R) return seg[id];
        int mid = (l + r) >> 1;
        int left = query(L, R, l, mid, id << 1);
        int right = query(L, R, mid + 1, r, id << 1 | 1);
        if (left == -1) return right;
        if (right == -1) return left;
        return h[left] < h[right] ? left : right;
    }
    int lca(int u, int v) {
        int l = first[u], r = first[v];
        if (l > r) swap(l, r);
        return query(l, r, 0, euler.size() - 1);
    }
};

struct LCA {
    int n, timer, LOG;
    vector<int> tin, tout;
    vector<vector<int>> up;

    LCA(const vector<vector<int>> &g, int root) {
        n = g.size(), timer = 0, LOG = ceil(log2(n));
        tin.resize(n), tout.resize(n);
        up.assign(n, vector<int>(LOG + 1));
        DFS(root, root, g);
    }
    void DFS(int u, int p, const vector<vector<int>> &g){
        tin[u] = ++timer, up[u][0] = p;
        for (int i = 1; i <= LOG; i++) up[u][i] = up[up[u][i - 1]][i - 1];
        for (int v : g[u]) if (v != p) DFS(v, u, g);
        tout[u] = ++timer;
    }
    bool is_ancestor(int u, int v) {
        return tin[u] <= tin[v] and tout[u] >= tout[v];
    }
    int lca(int u, int v) {
        if (is_ancestor(u, v)) return u;
        if (is_ancestor(v, u)) return v;
        for (int i = LOG; i >= 0; i--) 
            if (!is_ancestor(up[u][i], v)) u = up[u][i];
        return up[u][0];    
    }
};
