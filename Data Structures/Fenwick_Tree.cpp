#include <bits/stdc++.h>
using namespace std;
const int N = 1e6 + 10, M = 1e3 + 10;

/*--------------------------------------------------------------------------------------------------------
BIT :
Fenwick can only answer the queries of type [0, r] or [1, r] (prefix queries)
2 versions
0-base :
we use 0 base bit => bit[i] shows [g(i), i]
1-base :
we use 1 base bit => bit[i] shows (g(i), i]
Additional :
We can set all the bit to INF, for min queries
Range update Point query (just for sum queries) :
we make an arrey b, that if we get a partial sum of it, it becames a (it's not needed to do this, just for understanding the concept)
then for range update (l, r, x), we do point update (l, x) and (r + 1, -x)
and for point query (ind) we use answer (ind)
Range update Range query (just for sum queries) :
implemented (just for sum query)
NOTE that we always work with [ , ] segments in fenwick
--------------------------------------------------------------------------------------------------------*/

struct Data {
    int sum;
};

struct FenwickTree {
    int n;
    vector<Data> bit;

    void build(vector<int> &a) {
        n = a.size();
        bit.resize(n + 1);
        for (int i = 1; i <= n; i++) 
            bit[i].sum = 0;
        for (int i = 0; i < n; i++)
            update(i + 1, a[i]);
    }

    void update(int ind, int val) {
        for (; ind <= n; ind += ind & -ind) 
            bit[ind].sum += val;
    }
    
    int prefix_answer(int r) {
        int res = 0;
        for (; r > 0; r -= r & -r)   
            res += bit[r].sum;
        return res;    
    }

    int answer(int l, int r) {
        return prefix_answer(r) - prefix_answer(l - 1);
    }
};

struct FenwickTreeZeroBase {
    int n;
    vector<Data> bit;

    void build(vector<int> &a) {
        n = a.size();
        bit.resize(n);
        for (int i = 0; i < n; i++) 
            bit[i].sum = 0;
        for (int i = 0; i < n; i++)
            update(i, a[i]);
    }

    void update(int ind, int val) {
        for (; ind < n; ind = ind | (ind + 1)) 
            bit[ind].sum += val;
    }

    int prefix_answer(int r) {
        int res = 0;
        for (; r >= 0; r = (r & (r + 1)) - 1)   
            res += bit[r].sum;
        return res;    
    }
 
    int answer(int l, int r) {
        return prefix_answer(r) - prefix_answer(l - 1);
    }
};

struct FenwickTree2D {
    int n, m;
    vector<vector<Data>> bit;
 
    void build(vector<vector<int>> &a) {
        n = a.size(), m = a[0].size();
        bit.resize(n, vector<Data>(m));
        for (int i = 0; i < n; i++)
            for (int j = 0; j < m; j++)
                update(i, j, a[i][j]);
    }

    void update(int x, int y, int val) {
        for (int i = x; i < n; i = i | (i + 1)) 
            for (int j = y; j < m; j = j | (j + 1)) 
                bit[i][j].sum += val;
    }

    int prefix_answer(int x, int y) {
        int res = 0;
        for (int i = x; i >= 0; i = (i & (i + 1)) - 1)
            for (int j = y; j >= 0; j = (j & (j + 1)) - 1)
                res += bit[i][j].sum;
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

    void build(vector<int> &a) {
        n = a.size();
        vector<int> temp(n, 0);
        B1.build(temp);
        B2.build(temp);
        for (int i = 0; i < n; i++) 
            update(i + 1, i + 1, a[i]);
    }

    void update(int l, int r, int val) {
        B1.update(l, val);
        B1.update(r + 1, -val);
        B2.update(l, val * (l - 1));
        B2.update(r + 1, -val * r);
    }

    int answer(vector<Data> &bit, int r) {
        int res = 0;
        for (; r > 0; r -= r & -r) 
            res += bit[r].sum;
        return res;    
    }

    int prefix_answer(int r) {
        return answer(B1.bit, r) * r - answer(B2.bit, r);
    }

    int range_answer(int l, int r) {
        return prefix_answer(r) - prefix_answer(l - 1);
    }
};
