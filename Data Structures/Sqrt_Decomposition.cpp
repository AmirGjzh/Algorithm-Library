#include <bits/stdc++.h>
using namespace std;
using ll = long long int;
constexpr int N = 1e6 + 10, SQRT = 1e3 + 10;

/*============================================================================================================
Sqrt Decomposition (basic block-based range queries):

Idea:
  • Split the array into blocks of size B ≈ √n
  • Precompute an aggregate value per block
  • Answer range queries by:
    – Taking whole blocks when possible
    – Scanning elements only at the boundaries

What this implementation supports:
  • Point update: add a value to a single index
  • Range query: sum over [l, r] in O(√n)
  • Easily generalizable to any associative operation (min, max, gcd, xor, etc.)

Notes:
  • Indices are 0-based, ranges are inclusive
  • Block size is chosen as ceil(sqrt(n))
  • For range updates, lazy values per block can be added (not implemented here)

Complexity:
  • Build: O(n)
  • Point update: O(1)
  • Range query: O(√n)

Mo’s Algorithm (offline range queries):

Idea:
  • Reorder queries so that the range [L, R] changes slowly between queries
  • Maintain a sliding window and update a data structure as L/R move

Requirements:
  • Offline: all queries known in advance
  • No updates to the array
  • Must define:
    – add(index)
    – remove(index)
    – get_answer()

This implementation:
  • Answers queries based on element frequencies
  • Maintains a set ordered by (frequency, value)
  • Returns the value with maximum frequency in the current range

Complexity:
  • O((N + Q)·√N·cost(add/remove))

SqrtTree (advanced static range queries):

What it is:
  • A static range-query structure with O(1) query time
  • Supports any associative operation
  • Allows point updates in O(√n)

How it works (high level):
  • Recursively split the array into power-of-two blocks
  • Precompute prefix and suffix results inside blocks
  • Precompute results for spans crossing block boundaries
  • Query answers are assembled from at most 3 precomputed pieces

Notes:
  • Much more complex than Sparse Table or basic sqrt-decomposition
  • Useful when you need:
    – Very fast queries
    – Some level of updates
  • Overkill unless constraints demand it

Complexity:
  • Build: O(n·log(log(n)))
  • Query: O(1)
  • Update: O(√n)

Two Common √-Decomposition Patterns (theory):

1) Light / Heavy (√-Split)
  • Split updates by cost:
    – Light updates: apply immediately
    – Heavy updates: defer and handle during queries
  • Each query processes only O(n / √n) heavy updates

2) Batch Rebuild
  • Collect updates in batches of size B ≈ √n
  • After each batch, rebuild the entire structure
  • Queries handle only recent (≤ B) updates individually

Why these work:
  • Both balance total work to roughly O(Q·√n)
  • Choose based on:
    – Predictable per-query cost (Light/Heavy)
    – Occasional expensive rebuilds (Batch)
============================================================================================================*/

struct Data {
    ll sum = 0;
};

struct SqrtDecomposition {
    int n, B;
    vector<Data> array, block;

    explicit SqrtDecomposition(const vector<ll> &a) {
        n = static_cast<int>(a.size()), B = static_cast<int>(ceil(sqrt(n)));
        array.resize(n), block.resize((n + B - 1) / B);
        for (int i = 0; i < n; i++) array[i].sum = a[i];
        for (int i = 0; i < (n + B - 1) / B; i++) block[i].sum = 0;    
        for (int i = 0; i < n; i++) block[i / B].sum += a[i];
    }
    void update(const int ind, const ll val) {
        array[ind].sum += val, block[ind / B].sum += val;
    }
    [[nodiscard]] ll query(const int l, const int r) const {
        ll res = 0;
        for (int i = l; i <= r; ) 
            if (i % B == 0 and i + B - 1 <= r) res += block[i / B].sum, i += B;
            else res += array[i++].sum;
        return res;    
    }
};

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
        answer.assign(static_cast<int>(queries.size()), 0);
        this->array = array;
        this->queries = queries;
        sort(this->queries.begin(), this->queries.end());
    }
    void add(const int ind) {
        data_structure.erase({frequency[array[ind]], array[ind]});
        frequency[array[ind]]++;
        data_structure.insert({frequency[array[ind]], array[ind]});
    }
    void remove(const int ind) {
        data_structure.erase({frequency[array[ind]], array[ind]});
        frequency[array[ind]]--;
        data_structure.insert({frequency[array[ind]], array[ind]});
    }
    ll get_answer() const {
        return data_structure.rbegin()->second;
    }
    vector<ll> solve() {
        int cur_l = 0, cur_r = -1;
        for (const auto &[l, r, ind] : queries) {
            while (cur_l > l) add(--cur_l);
            while (cur_r < r) add(++cur_r);
            while (cur_l < l) remove(cur_l++);
            while (cur_r > r) remove(cur_r--);
            answer[ind] = get_answer();
        }
        return answer;
    }
};

struct SqrtTree {
    vector<Data> v;
    int n, lg, index_sz;
    vector<int> clz, layers, on_layer;
    vector<vector<Data>> pref, suff, between;

    static int log2up(const int n) {
        int res = 0;
        while (1 << res < n) res++;
        return res;
    }
    static Data op(const Data &a, const Data &b) {
        Data res{};
        res.sum = a.sum + b.sum;
        return res;
    }
    void build_block(const int layer, const int l, const int r) {
        pref[layer][l] = v[l];
        for (int i = l + 1; i < r; i++) pref[layer][i] = op(pref[layer][i - 1], v[i]);
        suff[layer][r - 1] = v[r - 1];
        for (int i = r - 2; i >= l; i--) suff[layer][i] = op(v[i], suff[layer][i + 1]);  
    }
    void build_between(const int layer, const int l_bound, const int r_bound, const int between_offs) {
        const int bszlog = (layers[layer] + 1) >> 1, bcntlog = layers[layer] >> 1, 
        bsz = 1 << bszlog, bcnt = (r_bound - l_bound + bsz - 1) >> bszlog;
        for (int i = 0; i < bcnt; i++) {
            Data ans;
            for (int j = i; j < bcnt; j++) {
                Data add = suff[layer][l_bound + (j << bszlog)];
                ans = i == j ? add : op(ans, add);
                between[layer - 1][between_offs + l_bound + (i << bcntlog) + j] = ans;
            }
        }
    }
    void build_between_zero() {
        const int bszlog = (lg + 1) >> 1;
        for (int i = 0; i < index_sz; i++) v[n + i] = suff[0][i << bszlog];
        build(1, n, n + index_sz, (1 << lg) - n);    
    }
    void update_between_zero(const int bid) {
        const int bszlog = (lg + 1) >> 1;
        v[n + bid] = suff[0][bid << bszlog];
        update(1, n, n + index_sz, (1 << lg) - n, n + bid);
    }
    void build(const int layer, const int l_bound, const int r_bound, const int between_offs) {
        if (layer >= static_cast<int>(layers.size())) return;
        const int bsz = 1 << ((layers[layer] + 1) >> 1);
        for (int l = l_bound; l < r_bound; l += bsz) {
            const int r = min(l + bsz, r_bound);
            build_block(layer, l, r), build(layer + 1, l, r, between_offs);
        }
        if (layer == 0) build_between_zero();
        else build_between(layer, l_bound, r_bound, between_offs);
    }
    void update(const int layer, const int l_bound, const int r_bound, const int between_offs, const int x) {
        if (layer >= static_cast<int>(layers.size())) return;
        const int bszlog = (layers[layer] + 1) >> 1, bsz = 1 << bszlog, block_id = (x - l_bound) >> bszlog, 
        l = l_bound + (block_id << bszlog), r = min(l + bsz, r_bound);
        build_block(layer, l, r);
        if (layer == 0) update_between_zero(block_id);
        else build_between(layer, l_bound, r_bound, between_offs);
        update(layer + 1, l, r, between_offs, x);
    }
    Data query(const int l, const int r, const int between_offs, const int base) {
        if (l == r) return v[l];
        if (l + 1 == r) return op(v[l], v[r]);
        const int layer = on_layer[clz[(l - base) ^ (r - base)]], bszlog = (layers[layer] + 1) >> 1, 
        bcntlog = layers[layer] >> 1, l_bound = (((l - base) >> layers[layer]) << layers[layer]) + base, 
        l_block = ((l - l_bound) >> bszlog) + 1, r_block = ((r - l_bound) >> bszlog) - 1;
        Data ans = suff[layer][l];
        if (l_block <= r_block) {
            const Data add = layer == 0 ? query(n + l_block, n + r_block, (1 << lg) - n, n) 
            :between[layer - 1][between_offs + l_bound + (l_block << bcntlog) + r_block];
            ans = op(ans, add);
        }
        ans = op(ans, pref[layer][r]);
        return ans;
    }
    Data query(const int l, const int r) {
        return query(l, r, 0, 0);
    }
    void update(const int x, const Data &item) {
        v[x] = item;
        update(0, 0, n, 0, x);
    }
    explicit SqrtTree(const vector<Data> &a) {
        n = static_cast<int>(a.size()), lg = log2up(n), v = a, clz.resize(1 << lg), on_layer.resize(lg + 1);
        clz[0] = 0;
        for (int i = 1; i < clz.size(); i++) clz[i] = clz[i >> 1] + 1;
        int tlg = lg;
        while (tlg > 1) on_layer[tlg] = static_cast<int>(layers.size()), layers.push_back(tlg), tlg = (tlg + 1) >> 1;    
        for (int i = lg - 1; i >= 0; i--) on_layer[i] = max(on_layer[i], on_layer[i + 1]);
        const int between_layers = max(0, static_cast<int>(layers.size()) - 1), bszlog = (lg + 1) >> 1, bsz = 1 << bszlog;
        index_sz = (n + bsz - 1) >> bszlog;
        v.resize(n + index_sz);
        pref.assign(layers.size(), vector<Data>(n + index_sz));
        suff.assign(layers.size(), vector<Data>(n + index_sz));
        between.assign(between_layers, vector<Data>((1 << lg) + bsz));
        build(0, 0, n, 0);
    }
};
