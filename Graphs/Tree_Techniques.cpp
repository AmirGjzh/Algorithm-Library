#include <bits/stdc++.h>
using namespace std;
using ll = long long int;

/*============================================================================================================
Path Query Resolution on Trees via Euler Tour and LCA Decomposition                                       

  1. Query Rewriting                                                                                        
    • For each original query (u,v), compute w = LCA(u,v)                                                
    • Replace path‑query (u,v) by two “half‑queries”: (v,w) and (u,w)                                      
                                                                                                          
  2. Euler Tour Preparation                                                                                  
    • Perform a full Euler tour of the tree, recording entry times for each vertex                       
    • As you traverse, build two lists of edges (or vertices) in tour‑order:                                
      – Forward list: edges in the direction of the DFS (parent→child)                                     
      – Reverse list: edges in the backtrack direction (child→parent)                                     
                                                                                                             
  3. Indexing & Data Structures                                                                              
    • Map each vertex x to its first‑occurrence index in the forward list                                 
    • Maintain a range–query structure (e.g. Fenwick/BIT or segment tree) supporting:                       
      – Add and Delete operations on intervals                                                             
      – Aggregate (sum/count/max/etc.) over a given index range                                           
                                                                                                             
  4. Answering a “Half‑Query” (v,w)                                                                           
    • Let i_v = first‑occurrence of v, i_w = first‑occurrence of w in the forward list                    
    • Extract the aggregate on [i_w, i_v] in the forward‑list DS                                           
    • Remove (delete) the same interval [i_w, i_v] from the reverse‑list DS to isolate the v→w path      
                                                                                                             
  5. Combining Results                                                                                       
    • After processing both (v,w) and (u,w), merge their two “half‑path” answers to reconstruct            
      the full (u→v) query result.                                                                          
                                                                                                             
Time Complexity: O((n + Q)·log(n)) for n vertices and Q path‑queries.                                     
============================================================================================================*/

/*============================================================================================================
Heavy-Light Decomposition (HLD)

Overview:
    Decompose a rooted tree into "heavy" and "light" paths to support efficient
    path and subtree operations using a segment tree with lazy propagation
Core Components:
  • Euler Tour Order (pos[]):
    - Assigns each node a DFS entry index, subtree of u → contiguous range [pos[u], pos[u] + sz[u] - 1]
  • Heavy Path Decomposition:
    - Each node’s heavy child maxes out subtree size, yielding O(log(n)) chain jumps per root-to-node path
  • Segment Tree:
    - Supports range-add and range-sum in O(log(n)) with lazy propagation

Supported Operations:
  1) update_path(u, v, val):   add val to all nodes on path u → v      O(log²(n))
  2) query_path(u, v):         sum over all nodes on path u → v        O(log²(n))
  3) update_subtree(u, val):   add val to entire subtree rooted at u   O(log(n))
  4) query_subtree(u):         sum over entire subtree rooted at u     O(log(n))
  5) lca(u, v):                lowest common ancestor of u and v       O(log(n))

Time & Memory Complexity:
  - Preprocessing: O(n) for DFS + decompose
  - Path ops:      O(log²(n)) per update/query
  - Subtree ops:   O(log(n)) per update/query
  - Memory:        O(n) for vectors, O(4*n) for segment tree arrays

Usage Notes:
  - Root defaults to 0, to use a different root, call DFS(root) and decompos(root,root)
  - For edge-based queries, map each edge to its deeper endpoint’s pos[]
============================================================================================================*/

struct SegmentTree {
    int n;
    vector<ll> st, lazy;

    SegmentTree(int n): n(n), st(n << 2, 0), lazy(n << 2, 0) {}
    void apply(int node, int l, int r, ll val) {
        st[node] += (r - l + 1) * val;
        lazy[node] += val;
    }
    void push(int node, int l, int r) {
        if (lazy[node] != 0) {
            int mid = (l + r) >> 1;
            apply(node << 1, l, mid, lazy[node]);
            apply(node << 1 | 1, mid + 1, r, lazy[node]);
            lazy[node] = 0;
        }
    }
    void update(int node, int l, int r, int ql, int qr, ll val) {
        if (qr < l or r < ql) return;
        if (ql <= l and r <= qr) {apply(node, l, r, val); return;}
        push(node, l, r);
        int mid = (l + r) >> 1;
        update(node << 1, l, mid, ql, qr, val);
        update(node << 1 | 1, mid + 1, r, ql, qr, val);
        st[node] = st[node << 1] + st[node << 1 | 1];
    }
    ll query(int node, int l, int r, int ql, int qr) {
        if (qr < l or r < ql) return 0;
        if (ql <= l and r <= qr) return st[node];
        push(node, l, r);
        int mid = (l + r) >> 1;
        return query(node << 1, l, mid, ql, qr) +
            query(node << 1 | 1, mid + 1, r, ql, qr);
    }
    void update(int l, int r, ll val) {update(1, 0, n - 1, l, r, val);}
    ll query(int l, int r) {return query(1, 0, n - 1, l, r);}
};

struct HeavyLightDecomposition {
    SegmentTree *seg;
    int n, root, cur_pos;
    vector<vector<int>> G;
    vector<int> par, depth, heavy, head, pos, sz;

    HeavyLightDecomposition(const vector<vector<int>> &G): n(G.size()), G(G) {
        seg = new SegmentTree(n);
        root = 0; par.assign(n, -1); depth.assign(n, 0); sz.assign(n, 1);
        heavy.assign(n, -1); head.assign(n, 0); pos.assign(n, 0); 
        cur_pos = 0; DFS(root); decompos(root, root);
    }
    void DFS(int u) {
        int max_sz = 0;
        for (int v : G[u]) if (v != par[u]) {
            par[v] = u, depth[v] = depth[u] + 1;
            DFS(v); sz[u] += sz[v]; 
            if (sz[v] > max_sz) max_sz = sz[v], heavy[u] = v;
        }
    }
    void decompos(int u, int h) {
        head[u] = h; pos[u] = cur_pos++;
        if (heavy[u] != -1) decompos(heavy[u], h);
        for (int v : G[u]) if (v != par[u] and v != heavy[u]) 
        decompos(v, v);
    }
    int lca(int u, int v) {
        while (head[u] != head[v]) 
            if (depth[head[u]] > depth[head[v]]) u = par[head[u]];
            else v = par[head[v]];
        return ((depth[u] < depth[v]) ? u : v);    
    }
    void update_path(int u, int v, ll val) {
        while (head[u] != head[v]) 
            if (depth[head[u]] > depth[head[v]]) {
                seg->update(pos[head[u]], pos[u], val);
                u = par[head[u]];
            } else {
                seg->update(pos[head[v]], pos[v], val);
                v = par[head[v]];
            }
        int l = min(pos[u], pos[v]), r = max(pos[u], pos[v]);
        seg->update(l, r, val);
    }
    void update_subtree(int u, ll val) {
        seg->update(pos[u], pos[u] + sz[u] - 1, val);
    }
    ll query_subtree(int u) {
        return seg->query(pos[u], pos[u] + sz[u] - 1);
    }
    ll query_path(int u, int v) {
        ll res = 0;
        while (head[u] != head[v]) 
            if (depth[head[u]] > depth[head[v]]) {
                res += seg->query(pos[head[u]], pos[u]);
                u = par[head[u]];
            } else {
                res += seg->query(pos[head[v]], pos[v]);
                v = par[head[v]];
            }
        int l = min(pos[u], pos[v]), r = max(pos[u], pos[v]);
        res += seg->query(l, r);
        return res;
    }
};

/*============================================================================================================
DSU-on-Tree ("small-to-large" merging) Technique
 
Overview:
    DSU-on-tree optimizes subtree queries (e.g., counting distinct values in each subtree) by always merging the
    smaller child's data-structure into the larger one, Each element moves at most O(log(n)) times, yielding
    O(n.log(n)) total complexity instead of quadratic
 
Two Variants:
  1) Global DS + Rollback
    • One global data-structure supporting insert, query, and undo (rollback) operations
    • Process light (small) children first, rollback, then process heavy (large) child without rollback, then reapply lights
    • Requires rollback support, constants are ~2× higher
 
  2) Per-Node DS via Pointer-Swapping
    • Each node has its own DS pointer, Recurse lights fully, then heavy, steal heavy's DS by pointer swap
    • Merge light DS into heavy DS (always small→large) with linear passes over light nodes
    • No rollback needed, lower constants
 
Example Application: Counting Distinct Colors in Each Subtree
    Given a tree with `n` nodes, each colored with an integer, compute for every node u the number of distinct
    colors in the subtree rooted at u
============================================================================================================*/

struct DsuOnTree {
    int n;
    vector<int> sz, ans;
    vector<vector<int>> G;

    DsuOnTree(const vector<vector<int>> &G): n(G.size()), G(G) {
        sz.assign(n, 1);
        ans.assign(n, -1);
    }

    vector<int> col;
    map<int, int> DS;

    void DFS_size(int u, int p) {
        for (int v : G[u]) if (v != p) {
            DFS_size(v, u);
            sz[u] += sz[v];
        }
    }
    void add(int u, int p, int big_child) {
        /* Add `u` to DS */
        DS[col[u]]++;
        for (int v : G[u]) if (v != p and v != big_child) add(v, u, big_child);
    }
    void remove(int u, int p) {
        /* Remove `u` from DS */
        DS[col[u]]--;
        if (DS[col[u]] == 0) DS.erase(col[u]);
        for (int v : G[u]) if (v != p) remove(v, u);
    }
    void DFS_solve(int u, int p, bool keep) {
        int max_size = 0, big_child = -1;
        for (int v : G[u]) if (v != p and sz[v] > max_size) 
            max_size = sz[v], big_child = v;
        for (int v : G[u]) if (v != p and v != big_child) DFS_solve(v, u, false);
        if (big_child != -1) DFS_solve(big_child, u, true);
        add(u, p, big_child);
        /* Ready to answer about `u` */
        if (!keep) remove(u, p);
    }
};

struct DsuOnTree {
    int n;
    vector<int> sz, ans;
    vector<vector<int>> G;

    DsuOnTree(const vector<vector<int>> &G): n(G.size()), G(G) {
        sz.assign(n, 1);
        ans.assign(n, -1);
    }

    vector<int> col;
    vector<map<int, int>*> DS;

    void DFS_size(int u, int p) {
        for (int v : G[u]) if (v != p) {
            DFS_size(v, u);
            sz[u] += sz[v];
        }
    }
    void DFS_solve(int u, int p) {
        int max_size = 0, big_child = -1;
        for (int v : G[u]) if (v != p and sz[v] > max_size) 
            max_size = sz[v], big_child = v;
        if (big_child != -1) DFS_solve(big_child, u), DS[u] = DS[big_child];
        else DS[u] = new map<int, int>();
        for (int v : G[u]) if (v != p and v != big_child) {
            DFS_solve(v, u);
            /* Merge small children to `u` */
            for (auto [c, cnt] : (*DS[v])) (*DS[u])[c] += cnt;
        }
        /* Add `u` to DS[u] */
        (*DS[u])[col[u]]++;
        /* Ready to answer about `u` */
    }
};

/*============================================================================================================
Centroid Decomposition on Trees

Overview:
    A divide-and-conquer technique that recursively splits the tree into “centroid” roots,
    yielding a decomposition tree of height O(log(n)), Useful for answering path or subtree
    queries (e.g. distances, closest marked node, etc.) by aggregating contributions
    up the centroid tree

Typical Preprocessing:
  – Build the centroid tree by calling `CentroidDecomposition cd(G)`
  – The array `centroid_parent` encodes the new tree’s edges

Extensions:
  • Distance arrays: to support distance-based queries, precompute for each centroid
    c a `dist[c][x] = distance(c → x)` for all x in c’s component (via one BFS/DFS)
  • Level array: optionally store `level[c]` = depth of c in the centroid tree
  • Data per centroid: e.g. a multiset or Fenwick tree at each c to maintain updates
    on its component
  • Complete update/answer: the stubs `update(int u)` and `answer(int u)` must be
    filled in with your application-specific logic (e.g. marking u “red” and querying
    closest red node)

Time & Memory Complexity:
  – Building:     O(n.log(n)) (each node participates in O(log(n)) levels of recursion)
  – Query/Update: O(log n) per call (walking centroid_root ← … ← root)
  – Memory:       O(n + sum of per-centroid data structures)

Usage Notes:
  – The centroid tree is rooted arbitrarily (here at the centroid of the original tree)
  – All “removed” flags are reset only during construction, after that, the original
    tree edges remain intact
  – Good for problems requiring dynamic path queries with offline or online updates
============================================================================================================*/

struct CentroidDecomposition {
    int n;
    vector<vector<int>> G;
    vector<bool> is_removed;
    vector<int> size, centroid_parent;

    CentroidDecomposition(const vector<vector<int>> &G): G(G) {
        n = (int) G.size();
        size.assign(n, 0);
        is_removed.assign(n, false);
        centroid_parent.assign(n, -1);
        decompos(0, -1);
    }
    int DFS_size(int u, int p) {
        size[u] = 1;
        for (int v : G[u]) if (v != p and !is_removed[v]) 
            size[u] += DFS_size(v, u);
        return size[u];
    }
    int centroid(int u, int p, int tree_size) {
        for (int v : G[u]) {
            if (v == p or is_removed[v]) continue;
            if (size[v] * 2 > tree_size) return centroid(v, u, tree_size);
        }
        return u;
    }
    void decompos(int u, int p) {
        int cent = centroid(u, -1, DFS_size(u, -1));
        /* solve the things that depends on `u`  */
        is_removed[cent] = true;
        centroid_parent[cent] = p;
        for (int v : G[cent]) if (!is_removed[v]) decompos(v, cent);
    }
    void update(int u) {
        /* propagate an “update” from node u up through all centroid ancestors */
        for (int v = u; v != -1; v = centroid_parent[v]) {
            
        }
    }
    int answer(int u) {
        /* combine contributions from all centroid ancestor */
        for (int v = u; v != -1; v = centroid_parent[v]) {
            
        }
    }
};
