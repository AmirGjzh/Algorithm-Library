#include <bits/stdc++.h>
using namespace std;
using ll = long long int;
const int N = 1e6 + 10, SQRT = 1e3 + 10;

/*============================================================================================================
Sqrt Decomposition (for sum queries, generalizable to any associative operation)

Description:
  • Partition an array of size n into ⌈n / B⌉ blocks of size B ≈ √n
  • Supports:
    – Point update + range sum query in O(√n) time
    – (With slight modifications) range‑update / point‑query and range‑update / range‑query
  • Generalizable: Replace “sum” with any associative combine (min, max, gcd, XOR, product, etc.) 
    by redefining how `array[i]` and `block[b]` combine values

Setup:
  • Let B = ⌈√n⌉, number of blocks = ⌈n / B⌉
  • Maintain:
    – `array[i]` = A[i] (individual elements)
    – `block[b]` = combine of A in block b

Usage:
  • `build(a)`  
    – Initialize `array` from input, build each block’s combined total in O(n)
  • `update(ind, v) 
    – O(1)
  • `query(l, r)`  
    – Combine over [l…r] by:
      • Taking whole blocks in O(√n)  
      • Scanning partial blocks at the ends (≤ 2B elements)

Notes:
  • All indices and segments are inclusive [l, r], 0-based
  • For **range‑update / point‑query**, store an extra `lazy[b]` per block:
    – On `update(l…r, x)`:  
      * Add `x` to `array[i]` for i in l…min(r, end_of_block(l/B))  
      * Add `x` to `lazy[b]` for each full block b in (l/B+1)…(r/B-1)  
      * Add `x` to `array[i]` for i in max(l, start_of_block(r/B))…r  
    – On point query(i): return `combine(array[i], lazy[i/B])`
  • For **range‑update / range‑query**, layer two sqrt‑decompositions or use BIT tricks

Time Complexity:
  • build: O(n)  
  • point update: O(1)  
  • range query: O(√n)  
  • range‑update / point‑query: O(√n) per update, O(1) per query  
  • range‑update / range‑query: O(√n) per operation with two decompositions
============================================================================================================*/

struct Data {
    ll sum;
};

struct SqrtDecomposition {
    int n, B;
    vector<Data> array, block;
    
    SqrtDecomposition(const vector<ll> &a) {
        n = int(a.size()), B = int(ceil(sqrt(n)));
        array.resize(n), block.resize((n + B - 1) / B);
        for (int i = 0; i < n; i++) array[i].sum = a[i];
        for (int i = 0; i < (n + B - 1) / B; i++) block[i].sum = 0;    
        for (int i = 0; i < n; i++) block[i / B].sum += a[i];
    }
    void update(int ind, ll val) {
        array[ind].sum += val, block[ind / B].sum += val;
    }
    ll query(int l, int r) {
        ll res = 0;
        for (int i = l; i <= r; ) 
            if (i % B == 0 and i + B - 1 <= r) res += block[i / B].sum, i += B;
            else res += array[i++].sum;
        return res;    
    }
};

/*============================================================================================================
Mo’s Algorithm

Description:
  • An offline technique to answer range‐query problems (no updates) in O((N + Q)·√N·F),
    where F is the cost of add/remove operations
  • Sorts queries by block of left index (size ≈ √N), then by right index

Key Steps:
 1. **Sort Queries** by `(l/S, r)` so that successive queries share most of their range
 2. **Maintain a Sliding Window** `[curL…curR]` and a custom DS:
    – add(idx)    adds `array[idx]` to your DS  
    – remove(idx) removes `array[idx]` from your DS  
 3. **Answer** each query from the DS (e.g. current mode)

Requirements:
  • Offline: all queries known in advance
  • DS must support add/remove in O(F)
  • Use inclusive [l,r] indexing (0‑based)
============================================================================================================*/

struct Query {
    int l, r, ind;
    bool operator<(const Query other) const {
        return make_pair(l / SQRT, r) < make_pair(other.l / SQRT, other.r);
    }
};

struct MoTechnique {
    vector<Query> queries;
    vector<ll> answer, array;
    set<pair<int, ll>> data_structure;
    unordered_map<ll, int> frequency;
    
    MoTechnique(const vector<Query> &queries, const vector<ll> &array) {
        answer.assign(int(queries.size()), 0);
        this->array = array;
        this->queries = queries;
        sort(this->queries.begin(), this->queries.end());
    }
    void add(int ind) {
        data_structure.erase({frequency[array[ind]], array[ind]});
        frequency[array[ind]]++;
        data_structure.insert({frequency[array[ind]], array[ind]});
    }
    void remove(int ind) {
        data_structure.erase({frequency[array[ind]], array[ind]});
        frequency[array[ind]]--;
        data_structure.insert({frequency[array[ind]], array[ind]});
    }
    ll get_answer() {
        return data_structure.rbegin()->second;
    }
    vector<ll> solve() {
        int cur_l = 0, cur_r = -1;
        for (const Query &q : queries) {
            while (cur_l > q.l) add(--cur_l);
            while (cur_r < q.r) add(++cur_r);
            while (cur_l < q.l) remove(cur_l++);
            while (cur_r > q.r) remove(cur_r--);
            answer[q.ind] = get_answer();
        }
        return answer;
    }
};

/*============================================================================================================
SqrtTree

Description:
  • A static (and partially dynamic) range‐query structure with:
    – O(n·log(log(n))) build time
    – O(1) query time for any associative op
    – O(√n) point‐update time
  • Works by recursively partitioning the array into power‑of‑two “layers”,
    precomputing prefix/suffix sums inside each block, plus “between” tables 
    for spans that cross block boundaries

Usage:
  • precomputation(a)  
    – Supply initial Data array `a`, builds all layers
  • query(l, r)  
    – Returns op(a[l], a[l+1], …, a[r]) in O(1)
  • update(x, newValue)  
    – Point‐update index x to `newValue` in O(√n), updating all affected blocks

Key Internals:
  – `layers[]`: which block‐size exponent live at each recursion depth  
  – `pref[h][i], suff[h][i]` = prefix/suffix combine in layer h’s blocks  
  – `between[h][i,j]` = combine of entire span covering block‐i to block‐j  
  – `clz[]` / `on_layer[]` = quick mapping from span‐length (log) to which layer to use  
============================================================================================================*/

struct SqrtTree {
    vector<Data> v;
    int n, lg, index_sz;
    vector<int> clz, layers, on_layer;
    vector<vector<Data>> pref, suff, between;

    inline int log2up(int n) {
        int res = 0;
        while ((1 << res) < n) res++;
        return res;
    }
    inline Data op(const Data &a, const Data &b) {
        Data res;
        res.sum = a.sum + b.sum;
        return res;
    }
    inline void build_block(int layer, int l, int r) {
        pref[layer][l] = v[l];
        for (int i = l + 1; i < r; i++) pref[layer][i] = op(pref[layer][i - 1], v[i]);
        suff[layer][r - 1] = v[r - 1];
        for (int i = r - 2; i >= l; i--) suff[layer][i] = op(v[i], suff[layer][i + 1]);  
    }
    inline void build_between(int layer, int l_bound, int r_bound, int between_offs) {
        int bszlog = (layers[layer] + 1) >> 1, bcntlog = layers[layer] >> 1, 
        bsz = 1 << bszlog, bcnt = (r_bound - l_bound + bsz - 1) >> bszlog;
        for (int i = 0; i < bcnt; i++) {
            Data ans;
            for (int j = i; j < bcnt; j++) {
                Data add = suff[layer][l_bound + (j << bszlog)];
                ans = (i == j) ? add : op(ans, add);
                between[layer - 1][between_offs + l_bound + (i << bcntlog) + j] = ans;
            }
        }
    }
    inline void build_between_zero() {
        int bszlog = (lg + 1) >> 1;
        for (int i = 0; i < index_sz; i++) v[n + i] = suff[0][i << bszlog];
        build(1, n, n + index_sz, (1 << lg) - n);    
    }
    inline void update_between_zero(int bid) {
        int bszlog = (lg + 1) >> 1;
        v[n + bid] = suff[0][bid << bszlog];
        update(1, n, n + index_sz, (1 << lg) - n, n + bid);
    }
    inline void build(int layer, int l_bound, int r_bound, int between_offs) {
        if (layer >= (int) layers.size()) return;
        int bsz = 1 << ((layers[layer] + 1) >> 1);
        for (int l = l_bound; l < r_bound; l += bsz) {
            int r = min(l + bsz, r_bound);
            build_block(layer, l, r), build(layer + 1, l, r, between_offs);
        }
        if (layer == 0) build_between_zero();
        else build_between(layer, l_bound, r_bound, between_offs);
    }
    inline void update(int layer, int l_bound, int r_bound, int between_offs, int x) {
        if (layer >= (int) layers.size()) return;
        int bszlog = (layers[layer] + 1) >> 1, bsz = 1 << bszlog, block_id = (x - l_bound) >> bszlog, 
        l = l_bound + (block_id << bszlog), r = min(l + bsz, r_bound);
        build_block(layer, l, r);
        if (layer == 0) update_between_zero(block_id);
        else build_between(layer, l_bound, r_bound, between_offs);
        update(layer + 1, l, r, between_offs, x);
    }
    inline Data query(int l, int r, int between_offs, int base) {
        if (l == r) return v[l];
        if (l + 1 == r) return op(v[l], v[r]);
        int layer = on_layer[clz[(l - base) ^ (r - base)]], bszlog = (layers[layer] + 1) >> 1, 
        bcntlog = layers[layer] >> 1, l_bound = (((l - base) >> layers[layer]) << layers[layer]) + base, 
        l_block = ((l - l_bound) >> bszlog) + 1, r_block = ((r - l_bound) >> bszlog) - 1;
        Data ans = suff[layer][l];
        if (l_block <= r_block) {
            Data add = (layer == 0) ? (query(n + l_block, n + r_block, (1 << lg) - n, n)) 
            :(between[layer - 1][between_offs + l_bound + (l_block << bcntlog) + r_block]);
            ans = op(ans, add);
        }
        ans = op(ans, pref[layer][r]);
        return ans;
    }
    inline Data query(int l, int r) {
        return query(l, r, 0, 0);
    }
    inline void update(int x, const Data &item) {
        v[x] = item;
        update(0, 0, n, 0, x);
    }
    SqrtTree(const vector<Data> &a) {
        n = a.size(), lg = log2up(n), v = a, clz.resize(1 << lg), on_layer.resize(lg + 1);
        clz[0] = 0;
        for (int i = 1; i < clz.size(); i++) clz[i] = clz[i >> 1] + 1;
        int tlg = lg;
        while (tlg > 1) on_layer[tlg] = layers.size(), layers.push_back(tlg), tlg = (tlg + 1) >> 1;    
        for (int i = lg - 1; i >= 0; i--) on_layer[i] = max(on_layer[i], on_layer[i + 1]);
        int between_layers = max(0, (int) layers.size() - 1), bszlog = (lg + 1) >> 1, bsz = 1 << bszlog;
        index_sz = (n + bsz - 1) >> bszlog;
        v.resize(n + index_sz);
        pref.assign(layers.size(), vector<Data>(n + index_sz));
        suff.assign(layers.size(), vector<Data>(n + index_sz));
        between.assign(between_layers, vector<Data>((1 << lg) + bsz));
        build(0, 0, n, 0);
    }
};

/*============================================================================================================
Two Sqrt Decomposition Patterns

1) Light/Heavy (√-Split) Technique 
Overview:
  • Classify each update as "light" if it costs ≤ B work, or "heavy" if it costs ≫ B
  • Apply light updates immediately, defer heavy updates by storing their parameter
  • On each query, incorporate all light effects (already applied) plus iterate over  
    at most O(n/B) heavy updates to compute their contributions

Key Steps:
  1. Choose threshold B ≈ √n (or √Q)
  2. For each update:
    – If cost ≤ B → perform update in O(cost)
    – Else → push parameters into heavy_list
  3. For each query:
    – Start with base state from applied light updates
    – Loop over heavy_list (size ≤ n/B) and apply each heavy effect in O(1) (or O(log(n)))
    – Return combined result

Why It Works:
  - Total work from light updates  = O(Q·B)
  - Per-query work from heavy list  = O(Q·(n/B))
  - Balanced at B ≈ √n → total ≈ O(Q.√n)
   
2) Batch-Rebuild (√-Blocks) Technique

Overview:
  • Group updates into batches of size B
  • After B updates, perform one global “rebuild” in O(n) to consolidate all effects
  • Between rebuilds, maintain a small list of “recent” updates and handle them individually on queries

Key Steps:
  1. Choose batch size B ≈ √m (or √n)
  2. Maintain:
    – main_state (reflects last rebuild)
    – recent_list of up to B new updates
  3. On each update:
    – Apply it logically to recent_list (O(1) or O(cost_recent))
    – If recent_list.size() == B → run rebuild():  
      • Consolidate all updates into main_state in O(n)  
      • Clear recent_list
  4. On each query:
    – Compute answer from main_state in O(1) (or O(log(n)))
    – For each update in recent_list (≤ B), adjust the answer in O(1) (or O(log(n)))

Why It Works:
  - Rebuilds cost O((m/B)·n)  
  - In-between queries cost O(m.B)  
  - Balanced at B ≈ √(n) → total ≈ O(n√m + m√m)

Comparison:
  – Light/Heavy: no big global steps, but each query pays for all deferred heavy events
  – Batch-Rebuild: periodic O(n) consolidations, queries only pay for small recent list  
  – Choose based on whether you prefer steady per-query cost vs. occasional big rebuilds
============================================================================================================*/
