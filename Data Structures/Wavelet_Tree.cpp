#include <bits/stdc++.h>
using namespace std;

/*--------------------------------------------------------------------------------------------------------
WaveletTree :
This is a data structure, like a segment tree(but on values not indices)
It can answer 3 types of queries :
1 - number of elments in [l, r] less than or equal to k
2 - number of occurrences of k in [l, r]
3 - kth smallest element in [l, r]
all the queries in O(Log(b))
In this implementation, we use array instead of vector
In build function, elements of array are in range [x, y]
Note that our query ranges are 1-base
Just this :)
--------------------------------------------------------------------------------------------------------*/

struct WaveletTree {
    int lo, hi;
    vector<int> b;
    WaveletTree *left, *right;

    void build(int *from, int *to, int x, int y) {
        lo = x, hi = y;
        if (lo == hi or from >= to) return;
        int mid = (lo + hi) >> 1;
        auto f = [mid](int x) {
            return x <= mid;
        };
        b.reserve(to - from + 1);
        b.push_back(0);
        for (auto it = from; it != to; it++) 
            b.push_back(b.back() + f(*it));
        auto pivot = stable_partition(from, to, f);
        left = new WaveletTree(), left->build(from, pivot, lo, mid);
        right = new WaveletTree(), right->build(pivot, to, mid + 1, hi);    
    }

    int kth(int l, int r, int k) {
        if (l > r) return 0;
        if (lo == hi) return lo;
        int in_left = b[r] - b[l - 1];
        int lb = b[l - 1], rb = b[r];
        if (k <= in_left) return this->left->kth(lb + 1, rb, k);
        return this->right->kth(l - lb, r - rb, k - in_left);
    }

    int LTE(int l, int r, int k) {
        if (l > r or k < lo) return 0;
        if (hi <= k) return r - l + 1;
        int lb = b[l - 1], rb = b[r];
        return this->left->LTE(lb + 1, rb, k) + this->right->LTE(l - lb, r - rb, k);
    }

    int count(int l, int r, int k) {
        if (l > r or k < lo or k > hi) return 0;
        if (lo == hi) return r - l + 1;
        int lb = b[l - 1], rb = b[r], mid = (lo + hi) >> 1;
        if (k <= mid) return this->left->count(lb + 1, rb, k);
        return this->right->count(l - lb, r - rb, k);
    }
};
