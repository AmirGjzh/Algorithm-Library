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

struct HeavyLight {
    SegmentTree *seg;
    int n, root, cur_pos;
    vector<vector<int>> G;
    vector<int> par, depth, heavy, head, pos, sz;

    HeavyLight(const vector<vector<int>> &G): n(G.size()), G(G) {
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
