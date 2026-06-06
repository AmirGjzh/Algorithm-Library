#include <bits/stdc++.h>
using namespace std;
using ll = long long int;

/*============================================================================================================
Fenwick Tree / Binary Indexed Tree (BIT)

Description:
  Fenwick Tree is a data structure that supports efficient prefix-based queries
  and point or range updates on an array using binary decomposition of indices.

Supported Operations:
  • Prefix sum query:
      sum(0..r) or sum(1..r)
  • Range sum query:
      sum(l..r) = prefix(r) − prefix(l−1)
  • Point update:
      A[i] += x

Indexing Variants:
  1. 1-based indexing (classic BIT)
     • bit[i] stores sum over range (i − LSB(i), i]
  2. 0-based indexing
     • bit[i] stores sum over range [i − LSB(i) + 1, i]
     • Uses bitwise transitions with (i | (i + 1))

Advanced Techniques:
  • Range Update + Point Query
    – Maintain a difference array D
    – A[i] = prefix_sum(D, i)
    – To add x to [l, r]:
        update(l, +x)
        update(r + 1, −x)

  • Range Update + Range Query (Two BIT Method)
    – Maintain two Fenwick trees B1 and B2
    – Prefix sum formula:
        sum(r) = prefix(B1, r) * r − prefix(B2, r)
    – Enables range add and range sum queries

Included Structures:
  1. FenwickTree
     • 1-based indexing
     • Point update + prefix / range sum

  2. FenwickTreeZeroBase
     • 0-based indexing variant
     • Same functionality with alternative index transitions

  3. FenwickTree2D
     • 2D Fenwick Tree
     • Point update and submatrix sum queries

  4. RangeUpdateRangeQuery
     • Uses two 1-based Fenwick trees
     • Supports range addition and range sum queries

Time Complexity:
  • 1D update/query: O(log(n))
  • 2D update/query: O(log(n)·log(m))
  • Build from array: O(n.log(n))

Memory Complexity:
  • 1D BIT: O(n)
  • 2D BIT: O(n·m)

Notes:
  • BIT supports operations forming an abelian group (e.g. sum)
  • Does not support min/max queries without modification
  • All range queries are inclusive: [l, r]
============================================================================================================*/

struct Data {
    ll sum;
};

struct FenwickTree {
    int n;
    vector<Data> bit;

    explicit FenwickTree(const vector<ll> &a) {
        n = static_cast<int>(a.size());
        bit.resize(n + 1, {});
        for (int i = 0; i < n; i++) update(i + 1, a[i]);
    }
    void update(int ind, const ll val) {
        for (; ind <= n; ind += ind & -ind) bit[ind].sum += val;
    }
    [[nodiscard]] ll prefix_answer(int r) const {
        ll res = 0;
        for (; r > 0; r -= r & -r) res += bit[r].sum;
        return res;    
    }
    [[nodiscard]] ll answer(const int l, const int r) const {
        return prefix_answer(r) - prefix_answer(l - 1);
    }
};

struct FenwickTreeZeroBase {
    int n;
    vector<Data> bit;

    explicit FenwickTreeZeroBase(const vector<ll> &a) {
        n = static_cast<int>(a.size());
        bit.resize(n, {});
        for (int i = 0; i < n; i++) update(i, a[i]);
    }
    void update(int ind, const ll val) {
        for (; ind < n; ind = ind | ind + 1) bit[ind].sum += val;
    }
    [[nodiscard]] ll prefix_answer(int r) const {
        ll res = 0;
        for (; r >= 0; r = (r & r + 1) - 1) res += bit[r].sum;
        return res;    
    }
    [[nodiscard]] ll answer(const int l, const int r) const {
        return prefix_answer(r) - prefix_answer(l - 1);
    }
};

struct FenwickTree2D {
    int n, m;
    vector<vector<Data>> bit;

    explicit FenwickTree2D(const vector<vector<ll>> &a) {
        n = static_cast<int>(a.size()), m = static_cast<int>(a[0].size());
        bit.resize(n, vector<Data>(m, {0}));
        for (int i = 0; i < n; i++) 
            for (int j = 0; j < m; j++) update(i, j, a[i][j]);
    }
    void update(const int x, const int y, const ll val) {
        for (int i = x; i < n; i = i | i + 1) 
            for (int j = y; j < m; j = j | j + 1) bit[i][j].sum += val;
    }
    [[nodiscard]] ll prefix_answer(const int x, const int y) const {
        ll res = 0;
        for (int i = x; i >= 0; i = (i & i + 1) - 1)
            for (int j = y; j >= 0; j = (j & j + 1) - 1) res += bit[i][j].sum;
        return res;
    }
    [[nodiscard]] ll answer(const int x1, const int y1, const int x2, const int y2) const {
        return prefix_answer(x2, y2) - prefix_answer(x1 - 1, y2) - 
        prefix_answer(x2, y1 - 1) + prefix_answer(x1 - 1, y1 - 1);
    }
};

struct RangeUpdateRangeQuery {
    int n;
    FenwickTree *B1, *B2;

    explicit RangeUpdateRangeQuery(const vector<ll> &a) {
        n = static_cast<int>(a.size());
        B1 = new FenwickTree(vector<ll>(n, 0));
        B2 = new FenwickTree(vector<ll>(n, 0));
        for (int i = 0; i < n; i++) update(i + 1, i + 1, a[i]);
    }
    void update(const int l, const int r, const ll val) const {
        B1->update(l, val), B1->update(r + 1, -val);
        B2->update(l, val * (l - 1)), B2->update(r + 1, -val * r);
    }
    [[nodiscard]] ll prefix_answer(const int r) const {
        return B1->prefix_answer(r) * r - B2->prefix_answer(r);
    }
    [[nodiscard]] ll range_answer(const int l, const int r) const {
        return prefix_answer(r) - prefix_answer(l - 1);
    }
};
