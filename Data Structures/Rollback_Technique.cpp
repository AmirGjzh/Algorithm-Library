#include <bits/stdc++.h>
using namespace std;

/*============================================================================================================
Dynamic Connectivity via Segment-Tree-over-Time + Rollbackable DSU

What this code does
  • Solves dynamic connectivity problems offline, where edges are added and removed over time.
  • We want an answer for every time instant (e.g. number of connected components).
  • Deletions are hard for DSU, so we avoid them by converting time into a segment tree.
  • Each time index [0..Q-1] is a leaf of the segment tree.

Main idea (time as a segment tree)
  • An edge that exists during time interval [L, R] is added to all segment-tree nodes
    that fully cover this interval (O(log(Q)) nodes).
  • We DFS the segment tree:
    – When entering a node, apply all its edge additions to DSU.
    – If we reach a leaf (single time instant), record the current answer.
    – When leaving the node, rollback the DSU to restore the previous state.
  • Because of rollback, each leaf sees the DSU exactly as it was at that time.

How add / delete queries are handled (offline)
  • While reading queries, we only track ADD times.
  • For each edge:
    – On ADD: remember its start time.
    – On DELETE at time t: create interval [start, t-1].
    – If an edge is never deleted, its interval is [start, Q-1].
  • These intervals are what we insert into the segment tree.

Core structures
  • DsuSave:
    – Stores exactly what is needed to undo one union operation.
  • DsuWithRollbacks:
    – DSU without path compression (important for rollback).
    – unite(u,v) merges two components and pushes a DsuSave if merge happened.
    – rollback() undoes exactly one successful unite().
    – comps always stores current number of connected components.
  • Query:
    – Represents a single operation stored in the segment tree (here: an edge u–v).
    – modified records whether unite() actually changed DSU state (so we know if rollback is needed).
  • QueryTree:
    – Segment tree over time [0..Q-1].
    – Each node stores a list of Queries active over that time interval.
    – DFS applies queries on entry and rolls them back on exit.

Why rollback DSU works here
  • Each unite() either:
    – does nothing (already connected), or
    – changes DSU state and pushes one save record.
  • During DFS, every successful unite is rolled back exactly once.
  • No path compression is used, because it would destroy historical structure.

Complexity
  • Each interval is added to O(log(Q)) nodes.
  • Each unite / rollback is O(log(N)) (union by rank).
  • Total complexity: O((number of intervals)·log(Q)·log(N)).
  • Memory: O(Q·log(Q) + N).

Notes & pitfalls
  • Time indices are inclusive: intervals are [L, R].
  • Deletion at time t means the edge is inactive starting from t, so interval ends at t-1.
  • Always canonicalize edges (min(u,v), max(u,v)) when mapping adds/deletes.
  • Rollback correctness depends on: one unite → at most one rollback.
  • This version avoids templates and heap-owned lambdas to prevent memory leaks.

Generalization theory (no code here)
  • This pattern works for any data structure that:
    – Can apply an operation,
    – Push undo information onto a stack,
    – Roll back exactly one operation.
  • Examples: bipartite-check DSU, rollback segment trees, dynamic MST variants.
  • To generalize, only the DS logic changes; the segment-tree-over-time idea stays the same.
============================================================================================================*/

struct DsuSave {
    int u, rnku, v, rnkv;
    DsuSave(const int u, const int rnku, const int v, const int rnkv): u(u), rnku(rnku), v(v), rnkv(rnkv) {}
};

struct DsuWithRollbacks {
    int comps;
    stack<DsuSave> op;
    vector<int> par, rnk;

    explicit DsuWithRollbacks(const int n) {
        par.resize(n), rnk.resize(n);
        for (int i = 0; i < n; i++) par[i] = i, rnk[i] = 0;
        comps = n;    
    }
    int find_set(const int u) {
        return u == par[u] ? u : find_set(par[u]);
    }
    bool unite(int u, int v) {
        u = find_set(u), v = find_set(v);
        if(u != v) {
            op.emplace(u, rnk[u], v, rnk[v]), comps--;
            if (rnk[u] > rnk[v]) swap(u, v);
            par[u] = v;
            if (rnk[v] == rnk[u]) rnk[v]++;
            return true;    
        }    
        return false;
    }
    void rollback() {
        if (op.empty()) return;
        const DsuSave x = op.top(); op.pop(); comps++;
        par[x.u] = x.u, rnk[x.u] = x.rnku;
        par[x.v] = x.v, rnk[x.v] = x.rnkv;
    }
};

struct Query {
    int u, v;
    bool modified{};
    Query(const int u, const int v): u(u), v(v) {}
};

struct QueryTree {
    int Q;
    DsuWithRollbacks *dsu;
    vector<vector<Query>> tree;
 
    QueryTree(const int n, const int Q) : Q(Q) {
        dsu = new DsuWithRollbacks(n);
        tree.resize((Q + 1) << 2);
    }
    void add_to_tree(Query &op, const int L, const int R, const int l, const int r, const int id = 1) {
        if (L > R) return;
        if (L == l and r == R) {
            tree[id].push_back(op);
            return;
        }
        const int mid = (l + r) >> 1;
        add_to_tree(op, L, min(R, mid), l, mid, id << 1);
        add_to_tree(op, max(L, mid + 1), R, mid + 1, r, id << 1 | 1);
    }
    void add_query(Query op, const int l, const int r) {
        add_to_tree(op, l, r, 0, Q - 1);
    }
    void dfs(vector<int> &ans, const int l, const int r, const int id = 1) {
        for (Query &op : tree[id]) op.modified = dsu->unite(op.u, op.v);
        if (l == r) ans[l] = dsu->comps;
        else {
            const int mid = (l + r) >> 1;
            dfs(ans, l, mid, id << 1);
            dfs(ans, mid + 1, r, id << 1 | 1);
        }        
        for (const Query &op : tree[id]) if (op.modified) dsu->rollback();
    }
    vector<int> solve() {
        vector<int> ans(Q);
        dfs(ans, 0, Q - 1);
        return ans;
    }
};
