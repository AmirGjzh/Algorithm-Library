#include <bits/stdc++.h>
using namespace std;
using ll = long long int;
#define lid (id << 1)
#define rid (id << 1 | 1)
#define mid ((l + r) >> 1)

/*============================================================================================================
Segment Trees – Complete Practical Reference

What a Segment Tree really is
  • A binary tree over an index range [0…n−1]
  • Each node represents a contiguous segment
  • Leaves represent single elements
  • Internal nodes store an aggregate value of their children using an associative operation
  • Query and update always traverse O(log(n)) nodes

This file contains MULTIPLE segment tree variants.
They are NOT interchangeable. Each one exists because a different class of problems demands it.

1) Classic Segment Tree (Point Update + Range Query):

Purpose:
  • Static or semi-dynamic array
  • Point updates
  • Range queries (sum / min / max / xor / gcd / etc.)

Key ideas:
  • Each node stores a summary of its segment
  • Merge children using "combine(left, right)"
  • Updates propagate from leaf → root

Complexity:
  • Build: O(n)
  • Update: O(log(n))
  • Query: O(log(n))
  • Memory: ~4n nodes

When to use:
  • Most standard problems
  • When updates affect only one index
  • When operations are associative

Typical problems:
  • Range sum queries
  • Max/min in subarray
  • Frequency accumulation
  • Prefix/suffix tricks

2) Merge Sort Tree (Segment Tree with Sorted Containers):

Purpose:
  • Range queries that require ordering
  • Example: smallest element ≥ x in [l,r], count of elements ≤ x

Key ideas:
  • Each node stores a sorted container (multiset or vector)
  • Built bottom-up by merging children containers
  • Query splits range → binary search in O(log(n)) nodes

Complexity:
  • Build: O(n.log²(n)) (multiset merge)
  • Query: O(log²(n))
  • Update: O(log²(n)) if supported

When to use:
  • Static or lightly dynamic data
  • Value-based range queries

Typical problems:
  • K-th number in a range
  • Count greater/smaller than x
  • Offline order statistics

3) Lazy Propagation Segment Tree (Range Update + Range Query):

Purpose:
  • Heavy range updates
  • Avoid updating all elements explicitly

Key ideas:
  • Each node stores a "lazy" value meaning:
      “This update applies to my children, but I haven't pushed it yet”
  • Push lazy values down only when needed
  • Keeps complexity logarithmic

Supported operations:
  • Range add + range sum
  • Range assign + range min/max (with care)

Complexity:
  • Update: O(log(n))
  • Query: O(log(n))
  • Memory: ~4n with lazy flags

When to use:
  • Interval updates
  • Large number of operations

Typical problems:
  • Range increments
  • Interval painting
  • Difference array on steroids

4) Iterative Segment Tree (Bottom-Up, Cache-Friendly):

Purpose:
  • Faster constant factors
  • Works well with lazy propagation if implemented carefully

Key ideas:
  • Tree stored in array [N .. 2N)
  • Parent = i/2
  • Leaves stored contiguously

Advantages:
  • Faster than recursive version
  • Easier memory layout
  • Good for time-critical tasks

When to use:
  • Tight time limits
  • Large constraints

5) 2D Segment Tree:

Purpose:
  • Queries on submatrices
  • Point updates in 2D grid

Key ideas:
  • Segment tree over rows
  • Each node contains another segment tree over columns

Complexity:
  • Build: O(n·m·log(n)·log(m))
  • Query: O(log(n)·log(m))
  • Memory: huge

When to use:
  • Small grids only
  • No better alternative exists

Typical problems:
  • Submatrix sum
  • 2D range queries

6) Persistent Segment Tree (Immutable Versions):

Purpose:
  • Access historical versions
  • Queries on past states

Key ideas:
  • Path copying: only O(log(n)) new nodes per update
  • Old versions remain untouched
  • Tree is functional / immutable

Complexity:
  • Update: O(log(n)) time & memory
  • Query: O(log(n))

When to use:
  • Offline queries
  • Time travel problems
  • K-th order statistics

Typical problems:
  • K-th smallest in subarray
  • Versioned queries
  • Offline prefix queries

7) Dynamic / Implicit Segment Tree:

Purpose:
  • Extremely large coordinate range (e.g. up to 1e18)
  • Sparse updates

Key ideas:
  • Nodes created only when needed
  • No full array allocation
  • Works like a binary trie over ranges

Complexity:
  • Update/query: O(log(C)), C = coordinate range
  • Memory proportional to number of touched nodes

When to use:
  • Sparse data
  • Huge index space
============================================================================================================*/

struct Node {
    bool mark = false;
    ll sum = 0, lazy = 0;
    multiset<ll> elements;
};

struct SegmentTree {
    int n;
    vector<ll> array;
    vector<Node> seg;

    explicit SegmentTree(const vector<ll> &array) {
        n = static_cast<int>(array.size());
        seg.resize(n << 2, Node());
        this->array = array;
        build_tree(0, n - 1);
    }
    static Node make_node(const ll val) {
        Node res;
        res.sum = val;
        return res;
    }
    static Node combine(const Node &l, const Node &r) {
        Node res;
        res.sum = l.sum + r.sum;
        return res;
    }
    void build_tree(const int l, const int r, const int id = 1) {
        if (l == r) seg[id] = make_node(array[l]);
        else {
            build_tree(l, mid, lid), build_tree(mid + 1, r, rid);
            seg[id] = combine(seg[lid], seg[rid]);
        }
    }
    Node answer(const int L, const int R, const int l, const int r, const int id = 1) {
        if (L > R) return make_node(0);
        if (L == l and r == R) return seg[id];
        return combine(answer(L, min(R, mid), l, mid, lid), answer(max(L, mid + 1), R, mid + 1, r, rid));
    }
    void update(const int pos, const ll val, const int l, const int r, const int id = 1) {
        if (l == r) seg[id] = make_node(array[l] += val);
        else {
            if (pos <= mid) update(pos, val, l, mid, lid);
            else update(pos, val, mid + 1, r, rid);
            seg[id] = combine(seg[lid], seg[rid]);
        }
    }
};

struct MergeSortTree {
    int n;
    vector<ll> array;
    vector<Node> seg;

    explicit MergeSortTree(const vector<ll> &array) {
        n = static_cast<int>(array.size());
        seg.resize(n << 2, Node());
        this->array = array;
        build_tree(0, n - 1);
    }
    static Node make_node(const ll val) {
        Node res;
        res.elements.insert(val);
        return res;
    }
    static Node combine(const Node &l, const Node &r) {
        Node res;
        res.elements.insert(l.elements.begin(), l.elements.end());
        res.elements.insert(r.elements.begin(), r.elements.end());
        return res;
    }
    void build_tree(const int l, const int r, const int id = 1) {
        if (l == r) seg[id] = make_node(array[l]);
        else {
            build_tree(l, mid, lid), build_tree(mid + 1, r, rid);
            seg[id] = combine(seg[lid], seg[rid]);
        }
    }
    ll answer(const int L, const int R, const ll val, const int l, const int r, const int id = 1) {
        if (L > R) return LLONG_MAX;
        if (L == l and r == R) {
            if (const auto pos = seg[id].elements.lower_bound(val); pos != seg[id].elements.end()) return *pos;
            return LLONG_MAX;
        }
        return min(answer(L, min(R, mid), val, l, mid, lid), answer(max(L, mid + 1), R, val, mid + 1, r, rid));
    }
    void update(const int pos, const ll val, const int l, const int r, const int id = 1) {
        seg[id].elements.erase(seg[id].elements.find(array[pos])), seg[id].elements.insert(array[pos] + val);
        if (l == r) array[pos] += val;
        else if (pos <= mid) update(pos, val, l, mid, lid);
        else update(pos, val, mid + 1, r, rid);
    }
};

struct LazyPropagation {
    int n;
    vector<ll> array;
    vector<Node> seg;

    explicit LazyPropagation(const vector<ll> &array) {
        n = static_cast<int>(array.size());
        seg.resize(n << 2, Node());
        this->array = array;
        build_tree(0, n - 1);
    }
    static Node make_node(const ll val) {
        Node res;
        res.sum = val;
        return res;
    }
    static Node combine(const Node &l, const Node &r) {
        Node res;
        res.sum = l.sum + r.sum;
        return res;
    }
    void push(const int id, const int l, const int r) {
        if (seg[id].mark) {
            seg[lid].sum += seg[id].lazy * (mid - l + 1), seg[lid].lazy += seg[id].lazy, seg[lid].mark = true;
            seg[rid].sum += seg[id].lazy * (r - mid), seg[rid].lazy += seg[id].lazy, seg[rid].mark = true;
            seg[id].lazy = 0, seg[id].mark = false;
        }
    }
    void build_tree(const int l, const int r, const int id = 1) {
        if (l == r) seg[id] = make_node(array[l]);
        else {
            build_tree(l, mid, lid), build_tree(mid + 1, r, rid);
            seg[id] = combine(seg[lid], seg[rid]);
        }
    }
    Node answer(const int L, const int R, const int l, const int r, const int id = 1) {
        if (L > R) return make_node(0);
        if (L == l and R == r) return seg[id];
        push(id, l, r);
        return combine(answer(L, min(R, mid), l, mid, lid), answer(max(L, mid + 1), R, mid + 1, r, rid));
    }
    void update(const int L, const int R, const ll val, const int l, const int r, const int id = 1) {
        if (L > R) return;
        if (L == l and R == r) seg[id].sum += val * (R - L + 1), seg[id].lazy += val, seg[id].mark = true;
        else {
            push(id, l, r);
            update(L, min(R, mid), val, l, mid, lid), update(max(L, mid + 1), R, val, mid + 1, r, rid);
            seg[id].sum = seg[lid].sum + seg[rid].sum;
        }
    }
};

struct SegmentIterative {
    int n, N, h;
    vector<ll> array;
    vector<Node> seg;

    explicit SegmentIterative(const vector<ll> &array) {
        n = static_cast<int>(array.size());
        this->array = array;
        N = 1; h = 0; while (N < n) N <<= 1, h++;
        seg.resize(N << 1);
        for (int id = 0; id < n; id++)
        seg[N + id] = make_node(array[id]);
        for (int id = N - 1; id > 0; id--) {
            seg[id] = make_node(0); recal(id);
        }
    }
    static Node make_node(const ll val) {
        Node res;
        res.sum = val;
        return res;
    }
    void apply(const int id, const ll val, const int len) {
        seg[id].sum += val * len;
        if (id < N) seg[id].lazy += val;
    }
    void push(const int node) {
        for (int s = h; s > 0; s--) {
            if (const int id = node >> s; seg[id].lazy != 0) {
                apply(lid, seg[id].lazy, 1 << (s - 1));
                apply(rid, seg[id].lazy, 1 << (s - 1));
                seg[id].lazy = 0;
            }
        }
    }
    void recal(const int id) {
        const int len = 1 << (h - (31 - __builtin_clz(id)));
        seg[id].sum = seg[lid].sum + seg[rid].sum
        + seg[id].lazy * len;
    }
    void pull(int id) {
        for (id >>= 1; id > 0; id >>= 1) recal(id);
    }
    ll answer(int l, int r) {
        if (l > r) return 0;
        l += N, r += N;
        push(l), push(r);
        ll res = 0; for (; l <= r; l >>= 1, r >>= 1) {
            if (l & 1) res += seg[l++].sum;
            if (!(r & 1)) res += seg[r--].sum;
        }
        return res;
    }
    void update(int l, int r, const ll val) {
        if (l > r) return;
        l += N, r += N;
        const int ll = l;
        const int rr = r;
        push(ll), push(rr);
        for (int len = 1; l <= r; l >>= 1, r >>= 1, len <<= 1) {
            if (l & 1) apply(l++, val, len);
            if (!(r & 1)) apply(r--, val, len);
        }
        pull(ll), pull(rr);
    }
};

struct SegmentTree2D {
    int n, m;
    vector<vector<ll>> array;
    vector<vector<Node>> seg;

    explicit SegmentTree2D(const vector<vector<ll>> &array) {
        n = static_cast<int>(array.size()), m = static_cast<int>(array[0].size());
        seg.resize(n << 2, vector<Node>(m << 2));
        this->array = array;
        build_x(0, n - 1);
    }
    static Node make_node(const ll val) {
        Node res;
        res.sum = val;
        return res;
    }
    static Node combine(const Node &l, const Node &r) {
        Node res;
        res.sum = l.sum + r.sum;
        return res;
    }
    void build_y(const int lx, const int rx, const int idx, const int l, const int r, const int id = 1) {
        if (l == r)
            if (lx == rx) seg[idx][id] = make_node(array[lx][l]);
            else seg[idx][id] = combine(seg[idx << 1][id], seg[idx << 1 | 1][id]);
        else {
            build_y(lx, rx, idx, l, mid, lid), build_y(lx, rx, idx, mid + 1, r, rid);
            seg[idx][id] = combine(seg[idx][lid], seg[idx][rid]);
        }
    }
    void build_x(const int l, const int r, const int id = 1) {
        if (l != r) build_x(l, mid, lid), build_x(mid + 1, r, rid);
        build_y(l, r, id, 0, m - 1);
    }
    Node answer_y(const int ly, const int ry, const int idx, const int l, const int r, const int id = 1) {
        if (ly > ry) return make_node(0);
        if (ly == l and ry == r) return seg[idx][id];
        return combine(answer_y(ly, min(ry, mid), idx, l, mid, lid), answer_y(max(ly, mid + 1), ry, idx, mid + 1, r, rid));
    }
    Node answer(const int lx, const int rx, const int ly, const int ry, const int l, const int r, const int id = 1) {
        if (lx > rx) return make_node(0);
        if (lx == l and rx == r) return answer_y(ly, ry, id, 0, m - 1);
        return combine(answer(lx, min(rx, mid), ly, ry, l, mid, lid), answer(max(lx, mid + 1), rx, ly, ry, mid + 1, r, rid));
    }
    void update_y(const int x, const int y, const ll val, const int lx, const int rx, const int idx, const int l, const int r, const int id = 1) {
        if (l == r)
            if (lx == rx) seg[idx][id] = make_node(array[x][y] += val);
            else seg[idx][id] = combine(seg[idx << 1][id], seg[idx << 1 | 1][id]);
        else {
            if (y <= mid) update_y(x, y, val, lx, rx, idx, l, mid, lid);
            else update_y(x, y, val, lx, rx, idx, mid + 1, r, rid);
            seg[idx][id] = combine(seg[idx][lid], seg[idx][rid]);
        }
    }
    void update(const int x, const int y, const ll val, const int l, const int r, const int id = 1) {
        if (l != r)
            if (x <= mid) update(x, y, val, l, mid, lid);
            else update(x, y, val, mid + 1, r, rid);
        update_y(x, y, val, l, r, id, 0, m - 1);
    }
};

struct Vertex {
    ll sum;
    Vertex *left, *right;

    explicit Vertex(const ll val) : sum(val), left(nullptr), right(nullptr) {}
    Vertex(Vertex *left, Vertex *right) : sum(0), left(left), right(right) {
        if (left) sum += left->sum;
        if (right) sum += right->sum;
    }
};

struct PersistentSegmentTree {
    int n;
    vector<ll> array;
    vector<Vertex*> versions;

    explicit PersistentSegmentTree(const vector<ll> &array) {
        n = static_cast<int>(array.size());
        this->array = array;
        versions.push_back(build_tree(0, n - 1));
    }
    Vertex* build_tree(const int l, const int r) {
        if (l == r) return new Vertex(array[l]);
        return new Vertex(build_tree(l, mid), build_tree(mid + 1, r));
    }
    [[nodiscard]] ll answer(const int L, const int R, const int version) const {
        return answer(L, R, 0, n - 1, versions[version]);
    }
    static ll answer(const int L, const int R, const int l, const int r, const Vertex *v) {
        if (L > R) return 0;
        if (L == l and R == r) return v->sum;
        return answer(L, min(mid, R), l, mid, v->left) + answer(max(mid + 1, L), R, mid + 1, r, v->right);
    }
    void update(const int pos, const ll val, const int version) {
        versions.push_back(update(pos, val, 0, n - 1, versions[version]));
    }
    Vertex* update(const int pos, const ll val, const int l, const int r, const Vertex *v) {
        if (l == r) return new Vertex(array[l] += val);
        if (pos <= mid) return new Vertex(update(pos, val, l, mid, v->left), v->right);
        return new Vertex(v->left, update(pos, val, mid + 1, r, v->right));
    }
};

struct DynamicVertex {
    ll sum = 0, left_bound, right_bound;
    DynamicVertex *left_child = nullptr, *right_child = nullptr;

    DynamicVertex(const ll left_bound, const ll right_bound) {
        this->left_bound = left_bound;
        this->right_bound = right_bound;
    }
    void extend() {
        if (!left_child and left_bound < right_bound) {
            const ll m = (left_bound + right_bound) >> 1;
            left_child = new DynamicVertex(left_bound, m);
            right_child = new DynamicVertex(m + 1, right_bound);
        }
    }
    void update(const ll pos, const ll val) {
        extend();
        if (left_child) {
            if (pos <= left_child->right_bound) left_child->update(pos, val);
            else right_child->update(pos, val);
            sum = left_child->sum + right_child->sum;
        }
        else sum += val;
    }
    ll answer(const ll lq, const ll rq) {
        if (lq <= left_bound and right_bound <= rq) return sum;
        if (max(left_bound, lq) > min(right_bound, rq)) return 0;
        extend();
        return left_child->answer(lq, rq) + right_child->answer(lq, rq);
    }
};
