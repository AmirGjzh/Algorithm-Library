#include <bits/stdc++.h>
using namespace std;
const int N = 1e6 + 10, SQRT = 1e3 + 10;

/*--------------------------------------------------------------------------------------------------------
SQRT :
We use this DS to partition our data into blocks of size Sqrt(n)
this gives us O(Sqrt(n)) for each query
this implementation is just for sum queries
we can modify it for other types of queries
NOTE that we always work with [ , ] segments 
Range update Point query (just for sum queries) :
this time, block just saves the changes on its own block, not the sum
for query (l, r, x), we add x to all tail elements (Sqrt(n)) and we add x to all blocks in this segment (Sqrt(n))
for query (ind), we answer array[ind] + block[ind / SQ]
Range update Range query (just for sum queries) :
here we use both 2 last cases(Range update point query and reverse)
we use 2 SqrtDecomposition
--------------------------------------------------------------------------------------------------------*/

struct Data {
    int sum;
};

struct SqrtDecomposition {
    int n;
    vector<Data> array, block;
    
    void build(vector<int> &a) {
        n = a.size();
        array.resize(n), block.resize(SQRT);
        for (int i = 0; i < n; i++) 
            array[i].sum = a[i];
        for (int i = 0; i < SQRT; i++) 
            block[i].sum = 0;    
        for (int i = 0; i < n; i++) 
            block[i / SQRT].sum += a[i];
    }

    void update(int ind, int val) {
        array[ind].sum += val;
        block[ind / SQRT].sum += val;
    }

    int answer(int l, int r) {
        int res = 0;
        for (int i = l; i <= r; ) 
            if (i % SQRT == 0 and i + SQRT - 1 < r) {
                res += block[i / SQRT].sum;
                i += SQRT;
            }
            else {
                res += array[i].sum;
                i++;
            }
        return res;    
    }
};

/*--------------------------------------------------------------------------------------------------------
Mo's Algorithm :
this algorithm is GOD :)
Note that we don't have any updates, just questions about segments
Note that this technique is OFFLINE
We sort the queries in a intresting way
and answer them one by one
first we have an empty segment, for each query we want to handle, we add or remove elements to our
segment, so we can answer the query
So we need to create a Data Structure base on our question and query types, 
so that we be able to add or remove from this DS, easily
in this implementation we handle query of type mode
mode = the element that has the most frequency
NOTE that we always work with [ , ] segments 
Order = (N + Q).Sqrt(N).F  (F is order of add and remove functions)
--------------------------------------------------------------------------------------------------------*/

struct Query {
    int l, r, ind;
    bool operator < (Query other) const {
        return make_pair(l / SQRT, r) < 
               make_pair(other.l / SQRT, other.r);
    }
};

struct MoTechnique {
    vector<int> answer, array;
    vector<Query> queries;

    set<pair<int, int>> data_structure;
    map<int, int> frequency;
    
    void precomputation(vector<Query> &queries, vector<int> &array) {
        answer.assign(queries.size(), 0);
        sort(queries.begin(), queries.end());
        this->queries = queries;
        this->array = array;
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

    int get_answer() {
        return data_structure.rbegin()->second;
    }

    vector<int> solve() {
        int cur_l = 0, cur_r = -1;
        for (Query q : queries) {
            while (cur_l > q.l) add(--cur_l);
            while (cur_r < q.r) add(++cur_r);
            while (cur_l < q.l) remove(cur_l++);
            while (cur_r > q.r) remove(cur_r--);
            answer[q.ind] = get_answer();
        }
        return answer;
    }
};

/*--------------------------------------------------------------------------------------------------------
SqrtTree :
A nice data structure that can answer range queries of type (a[l] op a[l + 1] op ... op a[r])
that op is an operation that is this true "x op (y op z) = (x op y) op z" very fast
It's really like segtree
Orders :
build => O(n.Log(Log(n)))
update => O(Sqrt(n))
query => O(1)
We also can do updates
both point update and range update (lazy propagation)
We just implement point update
--------------------------------------------------------------------------------------------------------*/

struct SqrtTree {
    int n, lg, index_sz;
    vector<Data> v;
    vector<int> clz, layers, on_layer;
    vector<vector<Data>> pref, suff, between;

    inline int log2up(int n) {
        int res = 0;
        while ((1 << res) < n) res++;
        return res;
    }

    inline Data op(Data &a, Data &b) {
        Data res;
        res.sum = a.sum + b.sum;
        return res;
    }

    inline void build_block(int layer, int l, int r) {
        pref[layer][l] = v[l];
        for (int i = l + 1; i < r; i++) 
            pref[layer][i] = op(pref[layer][i - 1], v[i]);
        suff[layer][r - 1] = v[r - 1];
        for (int i = r - 2; i >= l; i--)
            suff[layer][i] = op(v[i], suff[layer][i + 1]);  
    }

    inline void build_between(int layer, int l_bound, int r_bound, int between_offs) {
        int bszlog = (layers[layer] + 1) >> 1;
        int bcntlog = layers[layer] >> 1;
        int bsz = 1 << bszlog;
        int bcnt = (r_bound - l_bound + bsz - 1) >> bszlog;
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
        for (int i = 0; i < index_sz; i++) 
            v[n + i] = suff[0][i << bszlog];
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
            build_block(layer, l, r);
            build(layer + 1, l, r, between_offs);
        }
        if (layer == 0) build_between_zero();
        else build_between(layer, l_bound, r_bound, between_offs);
    }

    inline void update(int layer, int l_bound, int r_bound, int between_offs, int x) {
        if (layer >= (int) layers.size()) return;
        int bszlog = (layers[layer] + 1) >> 1;
        int bsz = 1 << bszlog;
        int block_id = (x - l_bound) >> bszlog;
        int l = l_bound + (block_id << bszlog);
        int r = min(l + bsz, r_bound);
        build_block(layer, l, r);
        if (layer == 0) update_between_zero(block_id);
        else build_between(layer, l_bound, r_bound, between_offs);
        update(layer + 1, l, r, between_offs, x);
    }

    inline Data query(int l, int r, int between_offs, int base) {
        if (l == r) return v[l];
        if (l + 1 == r) return op(v[l], v[r]);
        int layer = on_layer[clz[(l - base) ^ (r - base)]];
        int bszlog = (layers[layer] + 1) >> 1;
        int bcntlog = layers[layer] >> 1;
        int l_bound = (((l - base) >> layers[layer]) << layers[layer]) + base;
        int l_block = ((l - l_bound) >> bszlog) + 1;
        int r_block = ((r - l_bound) >> bszlog) - 1;
        Data ans = suff[layer][l];
        if (l_block <= r_block) {
            Data add = (layer == 0) ? (
                query(n + l_block, n + r_block, (1 << lg) - n, n)
            ) : (
                between[layer - 1][between_offs + l_bound + (l_block << bcntlog) + r_block]
            );
            ans = op(ans, add);
        }
        ans = op(ans, pref[layer][r]);
        return ans;
    }

    inline Data query(int l, int r) {
        return query(l, r, 0, 0);
    }

    inline void update(int x, Data &item) {
        v[x] = item;
        update(0, 0, n, 0, x);
    }

    inline void precomputation(vector<Data> &a) {
        n = a.size(), lg = log2up(n), v = a, clz.resize(1 << lg), on_layer.resize(lg + 1);
        clz[0] = 0;
        for (int i = 1; i < clz.size(); i++) 
            clz[i] = clz[i >> 1] + 1;
        int tlg = lg;
        while (tlg > 1) {
            on_layer[tlg] = layers.size();
            layers.push_back(tlg);
            tlg = (tlg + 1) >> 1;
        }    
        for (int i = lg - 1; i >= 0; i--) 
            on_layer[i] = max(on_layer[i], on_layer[i + 1]);
        int between_layers = max(0, (int) layers.size() - 1);
        int bszlog = (lg + 1) >> 1;
        int bsz = 1 << bszlog;
        index_sz = (n + bsz - 1) >> bszlog;
        v.resize(n + index_sz);
        pref.assign(layers.size(), vector<Data>(n + index_sz));
        suff.assign(layers.size(), vector<Data>(n + index_sz));
        between.assign(between_layers, vector<Data>((1 << lg) + bsz));
        build(0, 0, n, 0);
    }
};
