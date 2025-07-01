#include <bits/stdc++.h>
using namespace std;

/*============================================================================================================
Fenwick / Binary Indexed Tree (BIT)

Description:
  Efficient support for:
    • Prefix sums: sum of [0…r] or [1…r]
    • Range sums: sum of [l…r] via two prefix queries

Variants:
  • 1‑based BIT (standard): bit[i] covers range (i - LSB(i), i]
  • 0‑based BIT: bit[i] covers range [i - LSB(i)+1, i]; alternative indexing

Use Cases:
  • Point update + range sum query
  • Range Update + Point Query:
    – Build a difference array concept: maintaining array D where
      A[i] = prefix_sum(D, i)
    – To apply +x on [l, r], do:
      update(l, +x);
      update(r + 1, -x);
    – Then A[ind] = prefix_sum(D, ind)
  • Range update + range query (via two BITs, supporting sum)
    – Use formula:
      sum(r) = prefix(B1, r) * r - prefix(B2, r)

Structures included:
  1. `FenwickTree`: 1-based point update/ prefix & range sum
  2. `FenwickTreeZeroBase`: 0-based equivalent
  3. `FenwickTree2D`: 2D point update / submatrix sum
  4. `RangeUpdateRangeQuery`: supports range add + range sum with two BITs

Time Complexity:
  • Update/query: O(log(n)) for 1D, O(log(n).log(m)) for 2D
  • Build: O(n.log(n))

Notes:
  • BIT supports only sum (or other group ops with inverses)
  • No min/range min support unless special variant used
  • 2D BIT uses O(n·m) memory and same logic with double loops
  • Range updates/range queries via two BITs: add strategy using B1 and B2
  • All segments work on **inclusive [l, r]** ranges
============================================================================================================*/

struct Data {
    int sum;
};

struct FenwickTree {
    int n;
    vector<Data> bit;

    void build(const vector<int> &a) {
        n = a.size();
        bit.resize(n + 1);
        for (int i = 1; i <= n; i++) bit[i].sum = 0;
        for (int i = 0; i < n; i++) update(i + 1, a[i]);
    }
    void update(int ind, int val) {
        for (; ind <= n; ind += ind & -ind) bit[ind].sum += val;
    }
    int prefix_answer(int r) {
        int res = 0;
        for (; r > 0; r -= r & -r) res += bit[r].sum;
        return res;    
    }
    int answer(int l, int r) {
        return prefix_answer(r) - prefix_answer(l - 1);
    }
};

struct FenwickTreeZeroBase {
    int n;
    vector<Data> bit;

    void build(const vector<int> &a) {
        n = a.size();
        bit.resize(n);
        for (int i = 0; i < n; i++) bit[i].sum = 0;
        for (int i = 0; i < n; i++) update(i, a[i]);
    }
    void update(int ind, int val) {
        for (; ind < n; ind = ind | (ind + 1)) bit[ind].sum += val;
    }
    int prefix_answer(int r) {
        int res = 0;
        for (; r >= 0; r = (r & (r + 1)) - 1) res += bit[r].sum;
        return res;    
    }
    int answer(int l, int r) {
        return prefix_answer(r) - prefix_answer(l - 1);
    }
};

struct FenwickTree2D {
    int n, m;
    vector<vector<Data>> bit;
 
    void build(const vector<vector<int>> &a) {
        n = a.size(), m = a[0].size();
        bit.resize(n, vector<Data>(m));
        for (int i = 0; i < n; i++) 
            for (int j = 0; j < m; j++) update(i, j, a[i][j]);
    }
    void update(int x, int y, int val) {
        for (int i = x; i < n; i = i | (i + 1)) 
            for (int j = y; j < m; j = j | (j + 1)) bit[i][j].sum += val;
    }
    int prefix_answer(int x, int y) {
        int res = 0;
        for (int i = x; i >= 0; i = (i & (i + 1)) - 1)
            for (int j = y; j >= 0; j = (j & (j + 1)) - 1) res += bit[i][j].sum;
        return res;
    }
    int answer(int x1, int y1, int x2, int y2) {
        return prefix_answer(x1, y1) - prefix_answer(x1, y2 - 1) - 
        prefix_answer(x2 - 1, y1) + prefix_answer(x2 - 1, y2 - 1);
    }
};

struct RangeUpdateRangeQuery {
    int n;
    FenwickTree B1, B2;

    void build(const vector<int> &a) {
        n = a.size();
        B1.build(vector<int>(n, 0));
        B2.build(vector<int>(n, 0));
        for (int i = 0; i < n; i++) update(i + 1, i + 1, a[i]);
    }
    void update(int l, int r, int val) {
        B1.update(l, val), B1.update(r + 1, -val);
        B2.update(l, val * (l - 1)), B2.update(r + 1, -val * r);
    }
    int prefix_answer(int r) {
        return B1.prefix_answer(r) * r - B2.prefix_answer(r);
    }
    int range_answer(int l, int r) {
        return prefix_answer(r) - prefix_answer(l - 1);
    }
};
