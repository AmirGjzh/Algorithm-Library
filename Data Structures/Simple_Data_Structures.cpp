#include <bits/stdc++.h>
using namespace std;
using ll = long long int;

/*============================================================================================================
MinStack

Description:
  • A stack data structure that supports retrieving the current minimum element in O(1) time

Implementation Details:
  • Internally uses a `stack<pair<int,int>>` where each element is:
    – `.first`  = actual pushed value
    – `.second` = minimum value in the stack at the time of this push

Operations:
  • push(x)  
  • pop()
  • get_min()

Time Complexity:
  • All operations are O(1)

Use Cases:
  • When you need stack behavior plus fast access to the minimum so far
============================================================================================================*/

struct MinStack {
    stack<pair<ll, ll>> st;

    void push(ll a) {
        ll new_min = st.empty() ? a : min(a, st.top().second);
        st.push({a, new_min});
    }
    void pop() {
        st.pop();
    }
    ll top() {
        return st.top().first;
    }
    ll get_min() {
        return st.top().second;
    }
};

/*============================================================================================================
MinQueue

Description:
  • A queue data structure supporting retrieving the current minimum element in amortized O(1) time
  • Internally uses two “MinStacks” (each a stack of {value, current_min}) to manage enqueues and dequeues

How it works:
  – **s1** (in-stack) holds newly pushed elements, each entry’s second component is the minimum of s1 up to that point
  – **s2** (out-stack) holds elements for popping, when s2 is empty and a pop is requested, we pour all of s1 into s2:
    • While moving, we recompute the running minimum for s2
  – `get_min()` examines the minima at the tops of s1 and s2 (whichever is non-empty) and returns their smaller

Operations:
  • push(x)
  • pop()
  • get_min()  
  • front()

Time Complexity:
  • Amortized O(1) per operation (each element moves at most once from s1 to s2)
============================================================================================================*/

struct MinQueue {
    stack<pair<ll, ll>> s1, s2;

    void push(ll x) {
        ll new_min = s1.empty() ? x : min(x, s1.top().second);
        s1.push({x, new_min});
    }
    void pop() {
        if (s2.empty())
            while (!s1.empty()) {
                ll x = s1.top().first;
                s1.pop();
                ll new_min = s2.empty() ? x : min(x, s2.top().second);
                s2.push({x, new_min});
            }
        s2.pop();
    }
    ll get_min() {
        if (s1.empty() or s2.empty()) return s1.empty() ? s2.top().second : s1.top().second;
        else return min(s1.top().second, s2.top().second);    
    }
    ll front() {
        if (s2.empty()) 
            while (!s1.empty()) {
                ll x = s1.top().first;
                s1.pop();
                ll new_min = s2.empty() ? x : min(x, s2.top().second);
                s2.push({x, new_min});
            }
        return s2.top().first;    
    }
};

/*============================================================================================================
Sparse Table

Description:
  • A static-range query structure for immutable arrays
  • Preprocesses in O(n·log(n)) and answers queries over any associative operation:
    – Idempotent (e.g. min, gcd, bitwise AND): query in O(1)
    – General (e.g. sum): query in O(log(n)) via binary lifting

Key Points:
  1. Precomputation (build):
    – `st[k][i]` stores combine over interval `[i … i+2ᵏ−1]`
    – Build levels k=0…LOG–1 in nested loops: O(n·log(n))
  2. Querying:
    – Idempotent → O(1):  
      Let k = ⌊log₂(r−l+1)⌋. 
      `ans = combine(st[k][l], st[k][r−2ᵏ+1])`
    – General → O(log(n)):  
      Decompose `[l…r]` into powers of two from largest to smallest, combining those blocks

Usage:
  • `build(arr)` once
  • Then call:
    – `query_idem(l, r)` for min/RMQ/GCD in O(1)
    – `query_log(l, r)` for sum/any associative op in O(log(n))

Notes:
  • Works on **inclusive [l, r]** ranges (0-based)
  • `combine(a,a)` must equal `a` for idempotent O(1) queries
============================================================================================================*/

struct SparseTable {
    int n, LOG;
    vector<vector<ll>> st;

    ll combine(ll a, ll b) {
        return a + b;
    }
    SparseTable(const vector<ll> &arr) {
        n = int(arr.size()), LOG = 32 - __builtin_clz(n);
        st.assign(LOG, vector<ll>(n));
        copy(arr.begin(), arr.end(), st[0].begin());
        for (int i = 1; i < LOG; i++)
            for (int j = 0; j + (1 << i) <= int(arr.size()); j++) 
                st[i][j] = combine(st[i - 1][j], st[i - 1][j + (1 << (i - 1))]);
    }
    ll query_log(int l, int r) {
        ll ans = 0;
        for (int i = LOG - 1; i >= 0; i--) 
            if ((1 << i) <= r - l + 1) ans = combine(ans, st[i][l]), l += (1 << i);
        return ans;    
    }
    ll query_idem(int l, int r) {
        int span = r - l + 1, k = 31 - __builtin_clz(span);
        return combine(st[k][l], st[k][r - (1 << k) + 1]);
    }
};

/*============================================================================================================
Disjoint Sparse Table

Description:
  • A static-range query structure for immutable arrays, supporting any associative operation 
    (sum, product, min, gcd, etc.) in O(1) time per query
  • Conceptually similar to a Sparse Table but “disjoint” on power‑of‑two blocks
  • Requires that the array length be padded up to the next power of two, filling extra slots 
    with the operation’s identity element (e.g., 0 for sum, 1 for product)

Key Points:
  1. Padding  
    – Let n₀ = arr.size(), choose n = 2ᵏ ≥ n₀. 
    – Extend arr to length n by appending `identity` values
  2. Precomputation (build)  
    – Allocate `st[h][0…n-1]` for h = 0…k
    – `st[0][i] = arr[i]`
    – For each h ≥ 1, block size = 2ʰ, half = 2ʰ⁻¹:  
      • In each block [b…b+2ʰ−1], compute prefix-sums to the right from middle and to the left from middle+1 
      • That fills st[h][b…b+2ʰ−1] so that any query straddling the midpoint uses two precomputed halves
  3. Query  
    – If l == r, return `st[0][l]`
    – Else let d = l ⊕ r, and h = floor(log₂(d))
    – The segments [l…(end of its half)] and [(start of r’s half)…r] both lie wholly within a block of size 2ʰ⁺¹  
    – Return `op(st[h+1][l], st[h+1][r])`

Operations:
  • `build(arr)` — O(n·log(n)) preprocessing  
  • `answer(l, r)` — O(1) query

Notes:
  • Uses inclusive [l, r] indexing (0-based) 
  • Must initialize `identity` and the `op` function appropriately  
============================================================================================================*/

struct DisjointSparseTable {
    vector<vector<ll>> st;

    ll op(ll a, ll b) {
        return a + b;
    }
    DisjointSparseTable(const vector<ll> &arr) {
        int n = 1, identity = 0;
        while (n < int(arr.size())) n <<= 1;
        st.assign(__lg(n) + 1, vector<ll>(n, identity)), copy(arr.begin(), arr.end(), st[0].begin());
        for (int h = 1, range; (range = (1 << h)) <= n; h++) {
            int half = range >> 1;
            for (int i = half; i < n; i += range) {
                st[h][i - 1] = st[0][i - 1];
                for (int j = i - 2; j >= i - half; j--) st[h][j] = op(st[h][j + 1], st[0][j]);
                st[h][i] = st[0][i];
                for (int j= i + 1; j < i + half; j++) st[h][j] = op(st[h][j - 1], st[0][j]);   
            }
        }
    }
    ll answer(int l, int r) {
        if (l == r) return st[0][l];
        int h = 32 - __builtin_clz(l ^ r);
        return op(st[h][l], st[h][r]);
    }
};
