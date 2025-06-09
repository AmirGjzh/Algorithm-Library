#include <bits/stdc++.h>
using namespace std;
const int N = 1e6 + 10, LOG = 25;

/*--------------------------------------------------------------------------------------------------------
Minimum Stack :
we can get the min of the stack in O(1)
we use stack<pair<T, T>> and stack.top.second is always the min
--------------------------------------------------------------------------------------------------------*/

struct MinStack {
    stack<pair<int, int>> st;

    void push(int a) {
        int new_min = st.empty() ? a : min(a, st.top().second);
        st.push({a, new_min});
    }

    void pop() {
        st.pop();
    }

    int get_min() {
        return st.top().second;
    }
};

/*--------------------------------------------------------------------------------------------------------
Minimum Queue :
we can get the min of the queue in O(1)
we use 2 stacks, and the else is more like the minimum stack
--------------------------------------------------------------------------------------------------------*/

struct MinQueue {
    stack<pair<int, int>> s1, s2;

    void push(int x) {
        int new_min = s1.empty() ? x : min(x, s1.top().second);
        s1.push({x, new_min});
    }

    void pop() {
        if (s2.empty())
            while (!s1.empty()) {
                int x = s1.top().first;
                s1.pop();
                int new_min = s2.empty() ? x : min(x, s2.top().second);
                s2.push({x, new_min});
            }
        s2.pop();
    }

    int get_min() {
        if (s1.empty() or s2.empty())
            return s1.empty() ? s2.top().second : s1.top().second;
        else 
            return min(s1.top().second, s2.top().second);    
    }
};

/*--------------------------------------------------------------------------------------------------------
Sparse Table :
Sparse Table is a data structure that allows answering range queries. It can answer some range queries in
O(Log(n)) and O(1) for range-minimum-queries
NOTE that we always work with [ , ] segments in sparse table
It is usefull for static data
1 - Precomputation : uses dynamic programming to solve 2-power lenght ranges in O(n.Log(n))
2 - Answer : we split the range into ranges with lenght of 2-power, and we combine them in O(Log(n))
--------------------------------------------------------------------------------------------------------*/

struct SparseTable {
    int st[LOG][N];

    int combine(int a, int b) {
        return a + b;
    }

    void build(vector<int> &arr) {
        copy(arr.begin(), arr.end(), st[0]);
        for (int i = 1; i < LOG; i++)
            for (int j = 0; j + (1 << i) <= (int) arr.size(); j++) 
                st[i][j] = combine(st[i - 1][j], st[i - 1][j + (1 << (i - 1))]);
    }

    int answer(int l, int r) {
        int ans = 0;
        for (int i = LOG - 1; i >= 0; i--) 
            if ((1 << i) <= r - l + 1) {
                ans = combine(ans, st[i][l]);
                l += (1 << i);
            }
        return ans;    
    }
};

/*--------------------------------------------------------------------------------------------------------
Disjoint Sparse Table :
similar to normal sparse table, but we can answer queries in O(1)
we can do answer on ranges that has associative function (sum, product, mid, gcd, ...)
This DS is really like SqrtTree(without updates)
Note that we must make the size of our input, to a power of 2
So we add some identity elements (1 in product or 0 in sum)
Done :)
--------------------------------------------------------------------------------------------------------*/

struct DisjointSparseTable {
    vector<vector<int>> st;

    int op(int a, int b) {
        return a + b;
    }

    void build(vector<int> &arr) {
        int n = 1, identity = 0;
        while (n < (int) arr.size()) n <<= 1;
        st.assign(__lg(n) + 1, vector<int>(n, identity));
        copy(arr.begin(), arr.end(), st[0].begin());
        for (int h = 1, range; (range = (1 << h)) <= n; h++) {
            int half = range >> 1;
            for (int i = half; i < n; i += range) {
                st[h][i - 1] = st[0][i - 1];
                for (int j = i - 2; j >= i - half; j--)
                    st[h][j] = op(st[h][j + 1], st[0][j]);
                st[h][i] = st[0][i];
                for (int j= i + 1; j < i + half; j++) 
                    st[h][j] = op(st[h][j - 1], st[0][j]);   
            }
        }
    }

    int answer(int l, int r) {
        if (l == r) return st[0][l];
        int h = 32 - __builtin_clz(l ^ r);
        return op(st[h][l], st[h][r]);
    }
};
