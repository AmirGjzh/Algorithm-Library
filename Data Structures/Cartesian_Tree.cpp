#include <bits/stdc++.h>
using namespace std;
using ll = long long int;
constexpr int N = 1e6 + 10;

static thread_local std::mt19937_64 rng(static_cast<uint64_t>(chrono::steady_clock::now().time_since_epoch().count()));
inline int rand_priority() {
    return std::uniform_int_distribution(1, INT_MAX)(rng);
}

/*============================================================================================================
Implicit Treap & Cartesian Tree Toolkit

1. Implicit Treap (Cartesian Tree with Randomized Priority):

Description:
  • An implicit treap is a randomized balanced binary tree that represents a
      dynamic array without storing explicit keys.
  • The in-order traversal defines the sequence order (0-based indexing).
  • Heap property is maintained on a random priority to ensure balance.
  • Supports efficient range queries and range updates via lazy propagation.

Core Principles:
  • Implicit indexing: position = size(left subtree)
  • Expected balance due to random priorities
  • Split & merge are the fundamental operations
  • Lazy propagation enables fast range updates

Supported Operations:
  • split(root, L, R, k):
    - Splits the treap into:
        L = elements [0 … k−1]
        R = elements [k … n−1]
  • merge(root, L, R):
    - Merges two treaps where all indices in L precede R
  • insert(pos, value):
    - Insert value at index pos
  • remove(pos):
    - Remove element at index pos
  • range add update [L, R]:
    - Add a constant to all values in the interval
  • range sum query [L, R]:
    - Returns the sum over the interval
  • optional range reverse (flag supported)

Lazy Propagation:
  • lazy + mark:
    - Delays addition updates to entire subtrees
  • reverse flag:
    - Swaps children and propagates reversal lazily
  • push():
    - Applies pending lazy operations before descending
  • pull():
    - Recomputes subtree metadata after modifications

Metadata Maintained per Node:
  • size   : number of nodes in subtree
  • sum    : sum of values in subtree
  • val    : node’s own value
  • lazy   : pending additive update
  • flags  : mark (lazy active), reverse (subtree reversed)

Implementations Provided:
  1) Pointer-based treap (Node*):
    • Easier to reason about
    • Dynamic memory allocation
    • Useful for clarity and debugging

  2) Array-based treap (static memory):
    • Faster and cache-friendly
    • No dynamic allocation
    • Suitable for performance-critical problems
    • Maximum size fixed by constant N

Complexity:
  • split / merge:        O(log(n)) expected
  • insert / remove:      O(log(n)) expected
  • range update/query:   O(log(n)) expected
  • space:                O(n)

Applications:
  • Dynamic arrays / ropes
  • Range sum / range update problems
  • Sequence manipulation with reversals
  • Competitive programming data structures

Notes:
  • Always call push() before accessing children
  • Always call pull() after structural changes
  • Indices are 0-based
  • Random priority is essential for balance

2. Max-Heap Cartesian Tree Construction:

Description:
  • Builds a max-Cartesian tree from an array A (0-based).
  • Tree properties:
    – Heap-ordered by value (parent ≥ children)
    – In-order traversal yields the original array
  • Constructed in linear time using a monotonic stack.

Complexity:
  • Time:   O(n)
  • Space:  O(n)

Applications:
  • Range maximum query decomposition
  • Histogram and skyline problems
  • Longest path / visibility problems
  • Tree representations of array dominance

Notes:
  • Deterministic construction
  • No randomness involved
  • Stable with respect to equal elements
============================================================================================================*/

struct Node {
    ll val, sum, lazy; 
    int size, priority; 
    bool mark, reverse;
    Node *left, *right, *par;

    explicit Node(const ll val) {
        left = right = par = nullptr;
        lazy = 0, size = 1, priority = rand_priority();
        mark = reverse = false, this->val = sum = val;
    }
};

struct TreapPtr {
    Node *root;

    explicit TreapPtr(const vector<ll> &a): root(nullptr) {
        for (const ll i : a) merge(root, root, new Node(i));
    }
    static int size(const Node *root) {
        return root ? root->size : 0;
    }
    static ll sum(const Node *root) {
        return root ? root->sum : 0;
    }
    static void push(Node *root) {
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
            root->sum += static_cast<ll>(size(root)) * root->lazy, root->val += root->lazy, root->lazy = 0;
            root->mark = false;
        }
    }
    static void pull(Node *root) {
        if (!root) return;
        push(root->left), push(root->right);
        root->size = size(root->left) + size(root->right) + 1;
        root->sum = sum(root->left) + sum(root->right) + root->val;
    }
    static void split(Node *root, Node *&left, Node *&right, const int key) {
        if (!root) return void(left = right = nullptr);
        push(root);    
        if (size(root->left) < key) {
            split(root->right, root->right, right, key - size(root->left) - 1);
            if (root->right) root->right->par = root;
            left = root; left->par = nullptr;
        }
        else {
            split(root->left, left, root->left, key);
            if (root->left) root->left->par = root;
            right = root; right->par = nullptr;
        }
        pull(root);    
    }
    static void merge(Node *&root, Node *left, Node *right) {
        push(left), push(right);
        if (!left or !right) return void(root = left ? left : right);
        if (left->priority < right->priority) {
            merge(left->right, left->right, right);
            if (left->right) left->right->par = left;
            root = left; root->par = nullptr;
        }
        else {
            merge(right->left, left, right->left);
            if (right->left) right->left->par = right;
            root = right; root->par = nullptr;
        }
        pull(root);           
    }
    void insert(const int ind, const ll val) {
        Node *l, *r, *node = new Node(val);
        split(root, l, r, ind);   
        merge(l, l, node);
        merge(root, l, r);
    }
    void remove(const int ind) {
        Node *l, *r;
        split(root, l, r, ind);
        split(r, root, r, 1);
        merge(root, l, r);
    }
    void update(const int L, const int R, const ll val) {
        Node *a = nullptr, *b = nullptr, *c = nullptr;
        split(root, a, b, L);
        split(b, b, c, R - L + 1);
        b->lazy += val, b->mark = true;
        merge(root, a, b);
        merge(root, root, c);
    }
    ll answer(const int L, const int R) {
        Node *a, *b, *c;
        split(root, a, b, L);
        split(b, b, c, R - L + 1);
        const ll res = b->sum;
        merge(root, a, b);
        merge(root, root, c);
        return res;
    }
    static void output(Node *root) {
        if (!root) return;
        push(root);
        output(root->left);
        cout << root->val;
        output(root->right);
    }
};

struct TreapArr {
    bool mark[N]{}, reverse[N]{};
    ll val[N]{}, sum[N]{}, lazy[N]{};
    int root = 0, node_cnt = 0, size[N]{}, priority[N]{}, lc[N]{}, rc[N]{}, par[N]{};

    explicit TreapArr(const vector<ll> &a) {
        for (const ll value : a) {
            const int id = new_node(value);
            merge(root, root, id);
        }
    }
    int new_node(const ll value) {
        const int u = ++node_cnt;
        val[u] = sum[u] = value, lazy[u] = 0, size[u] = 1, priority[u] = rand_priority(); 
        mark[u] = reverse[u] = false, lc[u] = rc[u] = par[u] = 0;
        return u;
    }
    [[nodiscard]] int get_size(const int u) const {
        return u ? size[u] : 0;
    }
    [[nodiscard]] ll get_sum(const int u) const {
        return u ? sum[u] : 0;
    }
    void push(const int u) {
        if (!u) return;
        if (reverse[u]) {
            swap(lc[u], rc[u]);
            if (lc[u]) reverse[lc[u]] ^= true;
            if (rc[u]) reverse[rc[u]] ^= true;
            reverse[u] = false;
        }
        if (mark[u]) {
            val[u] += lazy[u];
            sum[u] += static_cast<ll>(get_size(u)) * lazy[u];
            if (lc[u]) lazy[lc[u]] += lazy[u], mark[lc[u]] = true;
            if (rc[u]) lazy[rc[u]] += lazy[u], mark[rc[u]] = true;
            lazy[u] = 0;
            mark[u] = false;
        }
    }
    void pull(const int u) {
        if (!u) return;
        push(lc[u]), push(rc[u]);
        size[u] = get_size(lc[u]) + get_size(rc[u]) + 1;
        sum[u] = get_sum(lc[u]) + get_sum(rc[u]) + val[u];
    }
    void split(const int u, int &left, int &right, const int key) {
        if (!u) return void(left = right = 0);
        push(u);
        if (get_size(lc[u]) < key) {
            split(rc[u], rc[u], right, key - get_size(lc[u]) - 1);
            if (rc[u]) par[rc[u]] = u;
            left = u; par[left] = 0;
        }
        else {
            split(lc[u], left, lc[u], key);
            if (lc[u]) par[lc[u]] = u;
            right = u; par[right] = 0;
        }
        pull(u);
    }
    void merge(int &u, const int left, const int right) {
        push(left), push(right);
        if (!left or !right) {u = left ? left : right; if (u) par[u] = 0; return;}
        if (priority[left] < priority[right]) {
            merge(rc[left], rc[left], right);
            if (rc[left]) par[rc[left]] = left;
            u = left; par[u] = 0;
        } else {
            merge(lc[right], left, lc[right]);
            if (lc[right]) par[lc[right]] = right;
            u = right; par[u] = 0;
        }
        pull(u);
    }
    void insert(const int pos, const ll value) {
        int l = 0, r = 0;
        split(root, l, r, pos);
        const int id = new_node(value);
        merge(l, l, id);
        merge(root, l, r);
    }
    void remove(const int pos) {
        int l = 0, m = 0, r = 0;
        split(root, l, m, pos);
        split(m, m, r, 1);
        merge(root, l, r);
    }
    void update(const int L, const int R, const ll value) {
        int a = 0, b = 0, c = 0;
        split(root, a, b, L);
        split(b, b, c, R - L + 1);
        if (b) lazy[b] += value, mark[b] = true;
        merge(root, a, b);
        merge(root, root, c);
    }
    ll answer(const int L, const int R) {
        int a = 0, b = 0, c = 0;
        split(root, a, b, L);
        split(b, b, c, R - L + 1);
        const ll res = get_sum(b);
        merge(root, a, b);
        merge(root, root, c);
        return res;
    }
    void output(const int u = 1) {
        if (!u) return;
        push(u);
        output(lc[u]);
        cout << val[u];
        output(rc[u]);
    }
};

struct CartesianTree {
    int n;
    vector<int> parent;
    vector<vector<int>> children;

    int build(const vector<ll> &A) {
        n = static_cast<int>(A.size()), parent.assign(n, -1), children.assign(n, {});
        stack<int> st;
        for (int i = 0; i < n; i++) {
            int last = -1;
            while (!st.empty() and A[st.top()] < A[i]) { last = st.top(); st.pop(); }
            if (!st.empty()) { parent[i] = st.top(); children[st.top()].push_back(i); }
            if (last != -1) { parent[last] = i; children[i].push_back(last); }
            st.push(i);
        }
        int root = -1;
        for (int i = 0; i < n; i++) if (parent[i] == -1) { root = i; break; }
        return root;
    }
};
