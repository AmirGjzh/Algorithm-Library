#include <bits/stdc++.h>
using namespace std;

/*============================================================================================================
Dynamic Connectivity via Segment‑Tree‑over‑Time + Rollback‑DSU

Problem Context:
  • We receive q time‑stamped operations (edge additions/removals) on a dynamic graph
  • We want, at each moment, the number of connected components
  • Removing edges in a DSU is hard => we use an offline technique with time segmentation + rollback

Technique Overview:
  1. Build a segment tree over the time axis [0…q−1]
  2. Each edge that is “alive” between times L and R is added to O(log(q)) segment‑tree nodes covering [L…R]
  3. Traverse this tree via DFS:
    – Upon entering a node, perform all its DSU 'unite' operations, recording them on a rollback stack
    – If at a leaf (query time), record current component count
    – Upon exiting the node, rollback those unites to restore DSU state
 4. DSUWithRollbacks supports easy rollback by tracking changes on a stack

Time & Memory Complexity:
  • Each edge addition spans O(log(q)) nodes → O(q.log(q)) total DSU ops
  • Each unite/find is O(log(n)) (path compression omitted to preserve rollback integrity)
  • Overall: O((q.log(q)).α(n)) with α(n) ≈ inverse Ackermann, memory O(n + q.log(q))

Key Components:
  – `DsuSave`: snapshot of DSU state (parents & ranks) before union
  – `DsuWithRollbacks`: DSU variant with `unite()`, `rollback()`, and `comps` counter
  – `Query`: edge + flag (`united` means actual merge occurred, for rollback logic)
  – `QueryTree`: segment tree over time that:
    • `add_query(q, L, R)`: schedules a union on [L…R]
    • DFS traversal calls unite at nodes, records comp count at leaves, and undoes union before returning
============================================================================================================*/

struct DsuSave {
    int u, rnku, v, rnkv;
    DsuSave() {}
    DsuSave(int u, int rnku, int v, int rnkv)
        : u(u), rnku(rnku), v(v), rnkv(rnkv) {}
};

struct DsuWithRollbacks {
    int comps;
    stack<DsuSave> op;
    vector<int> par, rnk;

    DsuWithRollbacks() {}
    DsuWithRollbacks(int n) {
        par.resize(n), rnk.resize(n);
        for (int i = 0; i < n; i++) 
            par[i] = i, rnk[i] = 0;
        comps = n;    
    }

    int find_set(int u) {
        return (u == par[u] ? u : find_set(par[u]));
    }

    bool unite(int u, int v) {
        u = find_set(u); 
        v = find_set(v);
        if(u != v) {
            op.push(DsuSave(u, rnk[u], v, rnk[v]));
            comps--;
            if (rnk[u] > rnk[v])
                swap(u, v);
            par[u] = v;
            if (rnk[v] == rnk[u])
                rnk[v]++;
            return true;    
        }    
        return false;
    }

    void rollback() {
        if (op.empty())
            return;
        DsuSave x = op.top();
        op.pop();
        comps++;
        par[x.u] = x.u, rnk[x.u] = x.rnku;
        par[x.v] = x.v, rnk[x.v] = x.rnkv;
    }
};

struct Query {
    int u, v;
    bool united;
    Query() {}
    Query(int u, int v) : u(u), v(v) {}
};

struct QueryTree {
    int Q;
    vector<vector<Query>> tree;
    DsuWithRollbacks dsu;

    QueryTree() {}
    QueryTree(int n, int Q) : Q(Q) {
        dsu = DsuWithRollbacks(n);
        tree.resize((Q + 1) << 2);
    }

    void add_to_tree(int L, int R, Query &q, int l, int r, int id = 1) {
        if (L > R)
            return;
        if (L == l and r == R) {
            tree[id].push_back(q);
            return;
        }    
        int mid = (l + r) >> 1;
        add_to_tree(L, min(R, mid), q, l, mid, id << 1);
        add_to_tree(max(L, mid + 1), R, q, mid + 1, r, id << 1 | 1);
    }

    void add_query(Query q, int l, int r) {
        add_to_tree(l, r, q, 0, Q - 1);
    }

    void dfs(vector<int> &ans, int l, int r, int id = 1) {
        for (Query &q : tree[id]) 
            q.united = dsu.unite(q.u, q.v);
        if (l == r)
            ans[l] = dsu.comps;
        else {
            int mid = (l + r) >> 1;
            dfs(ans, l, mid, id << 1);
            dfs(ans, mid + 1, r, id << 1 | 1);
        }        
        for (Query q : tree[id]) 
            if (q.united)
                dsu.rollback();
    }

    vector<int> solve() {
        vector<int> ans(Q);
        dfs(ans, 0, Q - 1);
        return ans;
    }
};
