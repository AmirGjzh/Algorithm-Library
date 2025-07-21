#include <bits/stdc++.h>
using namespace std;
using ll = long long int;
const int N = 1e6 + 10;

/*============================================================================================================
Implicit Treap (Cartesian Tree)

Description:
  • A randomized BST + heap ("treap") that supports array-like operations:
    – point insert/delete at any position
    – range queries (sum, min, max, etc.)
    – range updates (add, set, reverse, etc.)
  • Powered by split/merge + lazy propagation for efficiency

Core Concepts:
  – Keys are implicit: defined by in-order index (0-based), not stored explicitly
  – Random priority ensures expected O(log(n)) time per operation via heap property
  – Maintains subtree metadata (size, sum) and lazy flags (assign, reverse)

Operations:
  * split(root, left, right, key):
    – Divides treap into [0..key-1] and [key..n-1], updating size and pushing laziness
  * merge(root, left, right):
    – Combines two treaps into one; left < right by indices, pushes laziness before merging
  * insert(pos, value): split at pos → merge(left, new node) → merge with right
  * remove(pos): split at pos; split right by 1 → discard middle → merge rest
  * range update (add x to [L..R]): split into a, b, c → mark b's lazy flag → merge back
  * range query: split into a, b, c → answer b.sum → merge back

Lazy Handling (`push`):
  • Reverse flag: swaps children and propagates toggles downstream
  • Add-update flag: accumulates in `lazy`, marks node, later applied in `pull`

Complexity:
  • All operations (split, merge, insert, delete, update, query) run in expected O(log(n)) time

Use Cases:
  • Dynamic sequence manipulation (e.g., ropes, substring reversals, interval updates/queries)

Notes:
  • Works on implicit indices, keys = size(left subtree).
  • Always `push` before traversing `left`/`right` to ensure correct metadata
  • `pull` updates `size` and `sum` based on children after modifications
============================================================================================================*/

struct Node {
    ll val, sum, lazy; 
    int size, priority; 
    bool mark, reverse;
    Node *left, *right, *par;
    Node(int val) {
        left = right = par = nullptr;
        lazy = 0, size = 1, priority = rand();
        mark = reverse = false, this->val = sum = val;
    }
};

struct Treap {
    Node *root;

    Treap(const vector<ll> &a): root(nullptr) {
        for (ll i : a) merge(root, root, new Node(i));
    }
    int size(Node *root) {
        return root ? root->size : 0;
    }
    ll sum(Node *root) {
        return root ? root->sum : 0;
    }
    void push(Node *root) {
        if (!root) return;
        if (root->reverse) {
            root->reverse = false;
            swap(root->left, root->right);
            if (root->left) root->left->reverse ^= true;
            if (root->right) root->right->reverse ^= true;
        }
        if (root->mark) {
            if (root->left) root->left->lazy += root->lazy, root->left->mark = true;
            if (root->right) root->right->lazy += root->lazy, root->right->mark = true;
            root->sum += ll(size(root)) * root->lazy, root->val += root->lazy, root->lazy = 0;
            root->mark = false;
        }
    }
    void pull(Node *root) {
        if (!root) return;
        push(root->left), push(root->right);
        root->size = size(root->left) + size(root->right) + 1;
        root->sum = sum(root->left) + sum(root->right) + root->val;
    }
    void split(Node *root, Node *&left, Node *&right, int key) {
        if (!root) return void(left = right = nullptr);
        push(root);    
        if (size(root->left) < key) {
            if (root->right) root->right->par = nullptr;
            split(root->right, root->right, right, key - size(root->left) - 1);
            if (root->right) root->right->par = root;
            left = root;
        }
        else {
            if (root->left) root->left->par = nullptr;
            split(root->left, left, root->left, key);
            if (root->left) root->left->par = root;
            right = root;
        }
        pull(root);    
    }
    void merge(Node *&root, Node *left, Node *right) {
        push(left), push(right);
        if (!left or !right) return void(root = left ? left : right);
        if (left->priority < right->priority) {
            if (left->right) left->right->par = nullptr;
            merge(left->right, left->right, right);
            if (left->right) left->right->par = left;
            root = left;
        }
        else {
            if (right->left) right->left->par = nullptr;
            merge(right->left, left, right->left);
            if (right->left) right->left->par = right;
            root = right; 
        }
        pull(root);           
    }
    void insert(int ind, ll val) {
        Node *l, *r, *node = new Node(val);
        split(root, l, r, ind);   
        merge(l, l, node);
        merge(root, l, r);
    }
    void remove(int ind) {
        Node *l, *r;
        split(root, l, r, ind);
        split(r, root, r, 1);
        merge(root, l, r);
    }
    void update(int L, int R, ll val) {
        Node *a, *b, *c;
        split(root, a, b, L);
        split(b, b, c, R - L + 1);
        b->lazy += val, b->mark = true;
        merge(root, a, b);
        merge(root, root, c);
    }
    ll answer(int L, int R) {
        Node *a, *b, *c;
        split(root, a, b, L);
        split(b, b, c, R - L + 1);
        ll res = b->sum;
        merge(root, a, b);
        merge(root, root, c);
        return res;
    }
    void output(Node *root) {
        if (!root) return;
        push(root);
        output(root->left);
        cout << root->val;
        output(root->right);
    }
};

struct Treap {
    bool mark[N], reverse[N];
    ll val[N], sum[N], lazy[N];
    int root = 0, node_cnt = 0, size[N], priority[N], lc[N], rc[N], par[N];

    Treap(const vector<ll> &a) {
        for (ll value : a) {
            int id = new_node(value);
            merge(root, root, id);
        }
    }
    int new_node(ll value) {
        int u = ++node_cnt;
        val[u] = sum[u] = value, lazy[u] = 0, size[u] = 1, priority[u] = rand(); 
        mark[u] = reverse[u] = false, lc[u] = rc[u] = 0, par[u] = 0;
        return u;
    }
    int get_size(int u) {
        return u ? size[u] : 0;
    }
    ll get_sum(int u) {
        return u ? sum[u] : 0;
    }
    void push(int u) {
        if (!u) return;
        if (reverse[u]) {
            swap(lc[u], rc[u]);
            if (lc[u]) reverse[lc[u]] ^= true;
            if (rc[u]) reverse[rc[u]] ^= true;
            reverse[u] = false;
        }
        if (mark[u]) {
            val[u] += lazy[u];
            sum[u] += ll(get_size(u)) * lazy[u];
            if (lc[u]) lazy[lc[u]] += lazy[u], mark[lc[u]] = true;
            if (rc[u]) lazy[rc[u]] += lazy[u], mark[rc[u]] = true;
            lazy[u] = 0;
            mark[u] = false;
        }
    }
    void pull(int u) {
        if (!u) return;
        push(lc[u]), push(rc[u]);
        size[u] = get_size(lc[u]) + get_size(rc[u]) + 1;
        sum[u] = get_sum(lc[u]) + get_sum(rc[u]) + val[u];
    }
    void split(int u, int &left, int &right, int key) {
        if (!u) return void(left = right = 0);
        push(u);
        if (get_size(lc[u]) < key) {
            if (rc[u]) par[rc[u]] = 0;
            split(rc[u], rc[u], right, key - get_size(lc[u]) - 1);
            if (rc[u]) par[rc[u]] = u;
            left = u;
        }
        else {
            if (lc[u]) par[lc[u]] = 0;
            split(lc[u], left, lc[u], key);
            if (lc[u]) par[lc[u]] = u;
            right = u;
        }
        pull(u);
    }
    void merge(int &u, int left, int right) {
        push(left), push(right);
        if (!left or !right) return void(u = left ? left : right);
        if (priority[left] < priority[right]) {
            if (rc[u]) par[rc[u]] = 0;
            merge(rc[left], rc[left], right);
            if (rc[u]) par[rc[u]] = u;
            u = left;
        } else {
            if (lc[u]) par[lc[u]] = 0;
            merge(lc[right], left, lc[right]);
            if (lc[u]) par[lc[u]] = u;
            u = right;
        }
        pull(u);
    }
    void insert(int pos, ll value) {
        int l, r;
        split(root, l, r, pos);
        int id = new_node(value);
        merge(l, l, id);
        merge(root, l, r);
    }
    void remove(int pos) {
        int l, m, r;
        split(root, l, m, pos);
        split(m, m, r, 1);
        merge(root, l, r);
    }
    void update(int L, int R, ll value) {
        int a, b, c;
        split(root, a, b, L);
        split(b, b, c, R - L + 1);
        lazy[b] += value, mark[b] = true;
        merge(root, a, b);
        merge(root, root, c);
    }
    ll answer(int L, int R) {
        int a, b, c;
        split(root, a, b, L);
        split(b, b, c, R - L + 1);
        ll res = get_sum(b);
        merge(root, a, b);
        merge(root, root, c);
        return res;
    }
    void output(int u = 1) {
        if (!u) return;
        push(u);
        output(lc[u]);
        cout << val[u];
        output(rc[u]);
    }
};

/*============================================================================================================
Max‑Heap Cartesian Tree Construction

Description:
  • Builds the max‑Cartesian tree of an array A (0‑indexed)
    – Tree is heap‑ordered (parent ≥ children) and in‐order gives original sequence
  • Returns the root index, parent of each node, and children adjacency

Applications:
  • Fast visibility/“mountain gliding” routes (longest root‑to‑leaf path)
  • Range‑maximum query tree decompositions
  • Geometric skyline or histogram problem.

Notes:
  • Uses a single stack pass in O(n) time
  • For each new index i, pop smaller elements—they become children of i
  • Then push i onto the stack, The remaining top (if any) is i’s parent

Order: O(n) time, O(n) memory
============================================================================================================*/

struct CartesianTree {
    int n;
    vector<int> parent;
    vector<vector<int>> children;

    int build(const vector<ll> &A) {
        n = A.size(), parent.assign(n, -1), children.assign(n, {});
        stack<int> st;
        for (int i = 0; i < n; i++) {
            int last = -1;
            while (!st.empty() and A[st.top()] < A[i]) {last = st.top(); st.pop();}
            if (!st.empty()) {parent[i] = st.top(); children[st.top()].push_back(i);}
            if (last != -1) {parent[last] = i; children[i].push_back(last);}
            st.push(i);
        }
        int root = -1;
        for (int i = 0; i < n; i++) if (parent[i] == -1) {root = i; break;}
        return root;
    }
};
