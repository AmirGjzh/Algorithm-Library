#include <bits/stdc++.h>
using namespace std;
using ll = long long int;

/*============================================================================================================
Wavelet Tree

Description:
  • Value-based segment tree for static arrays
  • Answers range queries by value, not by index
  • All queries use 1-based indices for [l, r]

Supported Queries (O(log(value range))):
  • kth(l, r, k): k-th smallest element in subarray [l, r]
  • LTE(l, r, x): count of elements ≤ x in [l, r]
  • count(l, r, x): occurrences of value x in [l, r]

Structure & Invariants:
  • Each node represents a value interval [lo, hi]
  • mid = (lo + hi) / 2 splits values, not positions
  • b[i] = how many of the first i elements go to the left child (≤ mid)
  • Elements are stably partitioned so relative order is preserved

Construction:
  • Build from pointer range [from, to) and value range [minVal, maxVal]
  • Recursively partitions by value
  • Time: O(n·log(value range))
  • Memory: O(n·log(value range))

Notes:
  • Works on immutable arrays
  • Requires knowing global min/max value
  • Efficient alternative to segment trees for order-statistic queries
============================================================================================================*/

struct WaveletTree {
    ll lo, hi;
    vector<int> b;
    WaveletTree *left = nullptr, *right = nullptr;

    WaveletTree(ll *from, ll *to, const ll x, const ll y) {
        lo = x, hi = y;
        if (lo == hi or from >= to) return;
        ll mid = (lo + hi) >> 1;
        auto f = [mid](const ll p) {return p <= mid;};
        b.reserve(to - from + 1), b.push_back(0);
        for (auto it = from; it != to; it++) b.push_back(b.back() + f(*it));
        const auto pivot = stable_partition(from, to, f);
        left = new WaveletTree(from, pivot, lo, mid);
        right = new WaveletTree(pivot, to, mid + 1, hi);    
    }
    [[nodiscard]] ll kth(const int l, const int r, const int k) const {
        if (l > r) return -1;
        if (lo == hi) return lo;
        const int in_left = b[r] - b[l - 1], lb = b[l - 1], rb = b[r];
        if (k <= in_left) return left->kth(lb + 1, rb, k);
        return right->kth(l - lb, r - rb, k - in_left);
    }
    [[nodiscard]] int LTE(const int l, const int r, const ll k) const {
        if (l > r or k < lo) return 0;
        if (hi <= k) return r - l + 1;
        const int lb = b[l - 1], rb = b[r];
        return left->LTE(lb + 1, rb, k) + right->LTE(l - lb, r - rb, k);
    }
    [[nodiscard]] int count(const int l, const int r, const ll k) const {
        if (l > r or k < lo or k > hi) return 0;
        if (lo == hi) return r - l + 1;
        const ll mid = (lo + hi) >> 1;
        const int lb = b[l - 1], rb = b[r];
        if (k <= mid) return left->count(lb + 1, rb, k);
        return right->count(l - lb, r - rb, k);
    }
};
