#include <bits/stdc++.h>
using namespace std;

/*--------------------------------------------------------------------------------------------------------
This technique is used for this type of problems :
Adding and removing in a Data Structure, and answering some information about it
When removing from a DS is not that easy, we use this technique
This is OFFLINE
We make a Segment Tree over time
By time, I mean queries
If we have q queries, so we have q + 1 times (initial state and after each query)
Now we know that every element lives in DS for a segment of time between additi0n and deletion (segment of segment tree)
So we add this elemet (or query) in that segment of our segment tree (Log(q) segments)
So if we travel from root to the i'th leaf, we add elemets that exactly are in our DS in the time of i'th query
Ok! now we run a DFS on our segment tree, when we see a node, we add all the elements(or queries) that are 
in this node, and when we leave it, we undo our changes
That means we have to make our DS, able to roolback in time to undo changes
We can do this with stack of operations, to save information about a ADD operation(like DsuSave) 
so that when we want to undo, we take our state before the ADD operation correctly

Here we implement a DSU with ability to rollback
and solve the Dynamic Connectivity problem
--------------------------------------------------------------------------------------------------------*/

struct DsuSave {
    int u, rnku, v, rnkv;
    DsuSave() {}
    DsuSave(int u, int rnku, int v, int rnkv)
        : u(u), rnku(rnku), v(v), rnkv(rnkv) {}
};

struct DsuWithRollbacks {
    int comps;
    vector<int> par, rnk;
    stack<DsuSave> op;

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
