#include <bits/stdc++.h>
using namespace std;
using ll = long long int;

/*============================================================================================================
Data Structures: MinStack, MinQueue, Sparse Table, Disjoint Sparse Table

Overview:
  This file implements several classic data structures for fast minimum retrieval
  and static range queries. All structures assume 0-based indexing and inclusive [l, r] ranges.

1) MinStack:

  • Stack supporting O(1) retrieval of the current minimum
  • Each element stores {value, min_so_far}
  • Operations: push, pop, top, get_min
  • Time: O(1) per operation
  • Space: O(n)

2) MinQueue:

  • Queue supporting amortized O(1) minimum queries
  • Implemented using two MinStacks:
    – s1 for incoming elements
    – s2 for outgoing elements
  • Minimum is min(top(s1), top(s2))
  • Operations: push, pop, front, get_min
  • Time: amortized O(1)

3) Sparse Table:

  • Static range query structure for immutable arrays
  • Preprocessing: O(n.log(n))
  • Queries:
    – Idempotent ops (min, gcd): O(1)
    – General associative ops (sum): O(log(n))
  • st[k][i] stores result over range [i, i + 2^k − 1]

4) Disjoint Sparse Table:

  • Static range query structure with true O(1) queries
  • Supports any associative operation
  • Array is padded to next power of two using identity element
  • Query uses highest differing bit of l and r (l XOR r)
  • Preprocessing: O(n.log(n))

Typical Usage:
  • Sliding window minimum → MinQueue
  • Stack minimum tracking → MinStack
  • Static RMQ / GCD → Sparse Table
  • Static sum/product queries → Disjoint Sparse Table
============================================================================================================*/

struct MinStack {
    stack<pair<ll, ll>> st;

    void push(ll a) {
        ll new_min = st.empty() ? a : min(a, st.top().second);
        st.emplace(a, new_min);
    }
    void pop() { st.pop(); }
    ll top() { return st.top().first; }
    ll get_min() { return st.top().second; }
};

struct MinQueue {
    stack<pair<ll, ll>> s1, s2;

    void push(ll x) {
        ll new_min = s1.empty() ? x : min(x, s1.top().second);
        s1.emplace(x, new_min);
    }
    void pop() {
        if (s2.empty())
            while (!s1.empty()) {
                ll x = s1.top().first;
                s1.pop();
                ll new_min = s2.empty() ? x : min(x, s2.top().second);
                s2.emplace(x, new_min);
            }
        s2.pop();
    }
    ll get_min() {
        if (s1.empty() or s2.empty()) return s1.empty() ? s2.top().second : s1.top().second;
        return min(s1.top().second, s2.top().second);
    }
    ll front() {
        if (s2.empty())
            while (!s1.empty()) {
                ll x = s1.top().first;
                s1.pop();
                ll new_min = s2.empty() ? x : min(x, s2.top().second);
                s2.emplace(x, new_min);
            }
        return s2.top().first;
    }
};

struct SparseTable {
    int n, LOG;
    vector<vector<ll>> st;

    static ll combine(const ll a, const ll b) { return a + b; }
    explicit SparseTable(const vector<ll> &arr) {
        n = static_cast<int>(arr.size()), LOG = 32 - __builtin_clz(n);
        st.assign(LOG, vector<ll>(n));
        copy(arr.begin(), arr.end(), st[0].begin());
        for (int i = 1; i < LOG; i++)
            for (int j = 0; j + (1 << i) <= static_cast<int>(arr.size()); j++)
                st[i][j] = combine(st[i - 1][j], st[i - 1][j + (1 << (i - 1))]);
    }
    [[nodiscard]] ll query_log(int l, const int r) const {
        ll ans = 0;
        for (int i = LOG - 1; i >= 0; i--)
            if (1 << i <= r - l + 1) ans = combine(ans, st[i][l]), l += 1 << i;
        return ans;
    }
    [[nodiscard]] ll query_idem(const int l, const int r) const {
        const int span = r - l + 1, k = 31 - __builtin_clz(span);
        return combine(st[k][l], st[k][r - (1 << k) + 1]);
    }
};

struct DisjointSparseTable {
    vector<vector<ll>> st;

    static ll op(const ll a, const ll b) { return a + b; }
    explicit DisjointSparseTable(const vector<ll> &arr) {
        int n = 1;
        constexpr int identity = 0;
        while (n < static_cast<int>(arr.size())) n <<= 1;
        st.assign(__lg(n) + 1, vector<ll>(n, identity));
        copy(arr.begin(), arr.end(), st[0].begin());
        for (int h = 1, range; (range = (1 << h)) <= n; h++) {
            const int half = range >> 1;
            for (int i = half; i < n; i += range) {
                st[h][i - 1] = st[0][i - 1];
                for (int j = i - 2; j >= i - half; j--) st[h][j] = op(st[h][j + 1], st[0][j]);
                st[h][i] = st[0][i];
                for (int j= i + 1; j < i + half; j++) st[h][j] = op(st[h][j - 1], st[0][j]);
            }
        }
    }
    [[nodiscard]] ll answer(const int l, const int r) const {
        if (l == r) return st[0][l];
        const int h = 32 - __builtin_clz(l ^ r);
        return op(st[h][l], st[h][r]);
    }
};
