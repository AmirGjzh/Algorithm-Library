#include <bits/stdc++.h>
using namespace std;
using ll = long long int;

/*============================================================================================================
Wavelet Tree

Description:
  • A value-based segment tree enabling queries over values (not indices)
  • Supports (in O(log(value-range))):
    1. Count elements ≤ k in [l, r]
    2. Count occurrences of k in [l, r]
    3. Find k-th smallest element in [l, r]
  • Uses 1-based query indices

Structure:
  - Each node covers a value range [lo, hi]
  - `b[i]` = number of elements from `from[0..i-1]` that go to the left child
  - Partitions `from..to` by value around mid = (lo + hi)/2
  - Recursively builds left/right subtrees

Build:
  - Time: O(n·log(M)), where M = hi - lo
  - Must supply array and its value range [x, y]
  - Call `build(from, to, minVal, maxVal)`, passing:
    `from`: pointer to the start of the array (0‑based)
    `to`: pointer to one-past-the-end
    `minVal`: minimum element value in array
    `maxVal`: maximum element value in array
============================================================================================================*/

struct WaveletTree {
    ll lo, hi;
    vector<int> b;
    WaveletTree *left = nullptr, *right = nullptr;

    WaveletTree(ll *from, ll *to, ll x, ll y) {
        lo = x, hi = y;
        if (lo == hi or from >= to) return;
        ll mid = (lo + hi) >> 1;
        auto f = [mid](ll x) {return x <= mid;};
        b.reserve(to - from + 1), b.push_back(0);
        for (auto it = from; it != to; it++) b.push_back(b.back() + f(*it));
        auto pivot = stable_partition(from, to, f);
        left = new WaveletTree(from, pivot, lo, mid);
        right = new WaveletTree(pivot, to, mid + 1, hi);    
    }
    ll kth(int l, int r, int k) {
        if (l > r) return -1;
        if (lo == hi) return lo;
        int in_left = b[r] - b[l - 1], lb = b[l - 1], rb = b[r];
        if (k <= in_left) return left->kth(lb + 1, rb, k);
        return right->kth(l - lb, r - rb, k - in_left);
    }
    int LTE(int l, int r, ll k) {
        if (l > r or k < lo) return 0;
        if (hi <= k) return r - l + 1;
        int lb = b[l - 1], rb = b[r];
        return left->LTE(lb + 1, rb, k) + right->LTE(l - lb, r - rb, k);
    }
    int count(int l, int r, ll k) {
        if (l > r or k < lo or k > hi) return 0;
        if (lo == hi) return r - l + 1;
        ll mid = (lo + hi) >> 1;
        int lb = b[l - 1], rb = b[r];
        if (k <= mid) return left->count(lb + 1, rb, k);
        return right->count(l - lb, r - rb, k);
    }
};
