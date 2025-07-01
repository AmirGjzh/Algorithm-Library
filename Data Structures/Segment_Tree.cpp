#include <bits/stdc++.h>
using namespace std;
#define lid (id << 1)
#define rid (id << 1 | 1)
#define mid ((l + r) >> 1) 
const int N = 1e6 + 10, INF = 1e9 + 10;

/*============================================================================================================
Segment Tree

1. Basic Segment Tree (Point Update + Range Query)
  • Stores aggregate info (sum, min, max, etc.) over array segments in a binary-tree structure
  • Build: O(n), Query/Update: O(log(n))
  • Node represents [l..r], merging children via associative 'combine'
  • Memory: O(n) (≈4n in array form)

2. Merge Sort Tree (Range Query with Sorted Lists)
  • Each node holds sorted multiset of its segment
  • Query: e.g. find the minimum ≥ x in a range
  • Build: O(n.log²(n)), Query: O(log²(n)), Update also O(log²(n)) if supported

3. Lazy Propagation (Range Update + Range Query)
  • Adds 'lazy' markers for deferred updates over segments
  • Query/Update: O(log(n)) with additional memory per node
  • Supports range additions/assignments with range queries like sum/min/max

4. 2D Segment Tree (Submatrix Query/Point Update)
  • Tree of trees: segment tree over rows, each node has its own segment tree over columns
  • Build: O(n·m·log(n)·log(m)), Query/Update: O(log(n)·log(m))

5. Persistent Segment Tree (Versioned Trees)
  • Enables access to any historical version via path-copying
  • Per update: O(log(n)) time & memory, Query on any version: O(log(n))
  • Useful for rollback queries, Kth-order statistics, prefix histories

6. Dynamic / Implicit Segment Tree
  • Supports large coordinate ranges or sparse data without full array allocation
  • Nodes created dynamically during updates
  • Query/Update: O(log(C)), where C is coordinate range

Usage Notes:
  • Replace 'sum' with any associative operation (min, max, xor, gcd, etc.)
  • Adjust `combine()` and node data to your operation
  • Choose structure based on your application:
    – Merge Sort Tree for static queries like "count ≥ x"
    – Lazy Propagation for heavy range updates
    – Persistent Tree for querying past versions
    – Implicit Tree for sparsely indexed data
============================================================================================================*/

struct Node {
    int sum = 0, lazy = 0;
    bool mark = 0;
    multiset<int> elements;
};

struct SegmentTree {
    int n;
    vector<Node> seg;
    vector<int> array;

    void build(vector<int> &array) {
        n = array.size();
        seg.resize(n << 2);
        this->array = array;
        build_tree(0, n - 1);
    }
    Node make_node(int val) {
        Node res;
        res.sum = val;
        return res;
    }
    Node combine(Node l, Node r) {
        Node res;
        res.sum = l.sum + r.sum;
        return res;
    }
    void build_tree(int l, int r, int id = 1) {
        if (l == r) seg[id] = make_node(array[l]);
        else {
            build_tree(l, mid, lid), build_tree(mid + 1, r, rid);
            seg[id] = combine(seg[lid], seg[rid]);
        }    
    }
    Node answer(int L, int R, int l, int r, int id = 1) {
        if (L > R) return make_node(0);
        if (L == l and r == R) return seg[id];
        return combine(answer(L, min(R, mid), l, mid, lid), answer(max(L, mid + 1), R, mid + 1, r, rid));        
    }
    void update(int pos, int val, int l, int r, int id = 1) {
        if (l == r) seg[id] = make_node(array[l] += val);
        else {
            if (pos <= mid) update(pos, val, l, mid, lid);
            else update(pos, val, mid + 1, r, rid);
            seg[id] = combine(seg[lid], seg[rid]);    
        }    
    }
};

struct MergeSortTree {
    int n;
    vector<Node> seg;
    vector<int> array;

    void build(vector<int> &array) {
        n = array.size();
        seg.resize(n << 2);
        this->array = array;
        build_tree(0, n - 1);
    }
    Node make_node(int val) {
        Node res;
        res.elements.insert(val);
        return res;
    }
    Node combine(Node l, Node r) {
        Node res;
        res.elements.insert(l.elements.begin(), l.elements.end());
        res.elements.insert(r.elements.begin(), r.elements.end());
        return res;
    }
    void build_tree(int l, int r, int id = 1) {
        if (l == r) seg[id] = make_node(array[l]);
        else {
            build_tree(l, mid, lid), build_tree(mid + 1, r, rid);
            seg[id] = combine(seg[lid], seg[rid]);
        }    
    }
    int answer(int L, int R, int val, int l, int r, int id = 1) {
        if (L > R) return INF;
        if (L == l and r == R) {
            auto pos = seg[id].elements.lower_bound(val);
            if (pos != seg[id].elements.end()) return *pos;
            return INF;    
        }
        return min(answer(L, min(R, mid), val, l, mid, lid), answer(max(L, mid + 1), R, val, mid + 1, r, rid));        
    }
    void update(int pos, int val, int l, int r, int id = 1) {
        seg[id].elements.erase(seg[id].elements.find(array[pos])), seg[id].elements.insert(array[pos] + val);
        if (l == r) array[pos] += val;
        else if (pos <= mid) update(pos, val, l, mid, lid);
        else update(pos, val, mid + 1, r, rid);  
    }
};

struct LazyPropagation {
    int n;
    vector<Node> seg;
    vector<int> array;

    void build(vector<int> &array) {
        n = array.size();
        seg.resize(n << 2);
        this->array = array;
        build_tree(0, n - 1);
    }
    Node make_node(int val) {
        Node res;
        res.sum = val;
        return res;
    }
    Node combine(Node l, Node r) {
        Node res;
        res.sum = l.sum + r.sum;
        return res;
    }
    void push(int id, int l, int r) {
        if (seg[id].mark) {
            seg[lid].sum += seg[id].lazy * (mid - l + 1), seg[lid].lazy += seg[id].lazy, seg[lid].mark = 1;
            seg[rid].sum += seg[id].lazy * (r - mid), seg[rid].lazy += seg[id].lazy, seg[rid].mark = 1;
            seg[id].lazy = 0, seg[id].mark = 0;
        }
    }
    void build_tree(int l, int r, int id = 1) {
        if (l == r) seg[id] = make_node(array[l]);
        else {
            build_tree(l, mid, lid), build_tree(mid + 1, r, rid);
            seg[id] = combine(seg[lid], seg[rid]);
        }    
    }
    Node answer(int L, int R, int l, int r, int id = 1) {
        if (L > R) return make_node(0);
        if (L == l and R == r) return seg[id];
        push(id, l, r);
        return combine(answer(L, min(R, mid), l, mid, lid), answer(max(L, mid + 1), R, mid + 1, r, rid));
    }
    void update(int L, int R, int val, int l, int r, int id = 1) {
        if (L > R) return;
        if (L == l and R == r) seg[id].sum += val * (R - L + 1), seg[id].lazy += val, seg[id].mark = 1;
        else {
            push(id, l, r);
            update(L, min(R, mid), val, l, mid, lid), update(max(L, mid + 1), R, val, mid + 1, r, rid);
            seg[id].sum = seg[lid].sum + seg[rid].sum;
        }
    }
};

struct SegmentTree2D {
    int n, m;
    vector<vector<Node>> seg;
    vector<vector<int>> array;

    void build(vector<vector<int>> &array) {
        n = array.size(), m = array[0].size();
        seg.resize(n << 2, vector<Node>(m << 2));
        this->array = array;
        build_x(0, n - 1);
    }
    Node make_node(int val) {
        Node res;
        res.sum = val;
        return res;
    }
    Node combine(Node l, Node r) {
        Node res;
        res.sum = l.sum + r.sum;
        return res;
    }
    void build_y(int lx, int rx, int idx, int l, int r, int id = 1) {
        if (l == r) 
            if (lx == rx) seg[idx][id] = make_node(array[lx][l]);
            else seg[idx][id] = combine(seg[idx << 1][id], seg[idx << 1 | 1][id]);    
        else {
            build_y(lx, rx, idx, l, mid, lid), build_y(lx, rx, idx, mid + 1, r, rid);
            seg[idx][id] = combine(seg[idx][lid], seg[idx][rid]);
        }
    } 
    void build_x(int l, int r, int id = 1) {
        if (l != r) build_x(l, mid, lid), build_x(mid + 1, r, rid);
        build_y(l, r, id, 0, m - 1);
    }
    Node answer_y(int ly, int ry, int idx, int l, int r, int id = 1) {
        if (ly > ry) return make_node(0);
        if (ly == l and ry == r) return seg[idx][id];
        return combine(answer_y(ly, min(ry, mid), idx, l, mid, lid), answer_y(max(ly, mid + 1), ry, idx, mid + 1, r, rid));
    }
    Node answer(int lx, int rx, int ly, int ry, int l, int r, int id = 1) {
        if (lx > rx) return make_node(0);
        if (lx == l and rx == r) return answer_y(ly, ry, id, 0, m - 1);
        return combine(answer(lx, min(rx, mid), ly, ry, l, mid, lid), answer(max(lx, mid + 1), rx, ly, ry, mid + 1, r, rid));    
    }
    void update_y(int x, int y, int val, int lx, int rx, int idx, int l, int r, int id = 1) {
        if (l == r) 
            if (lx == rx) seg[idx][id] = make_node(array[x][y] += val);
            else seg[idx][id] = combine(seg[idx << 1][id], seg[idx << 1 | 1][id]);
        else {
            if (y <= mid) update_y(x, y, val, lx, rx, idx, l, mid, lid);
            else update_y(x, y, val, lx, rx, idx, mid + 1, r, rid);
            seg[idx][id] = combine(seg[idx][lid], seg[idx][rid]);
        }
    }
    void update(int x, int y, int val, int l, int r, int id = 1) {
        if (l != r) 
            if (x <= mid) update(x, y, val, l, mid, lid);
            else update(x, y, val, mid + 1, r, rid);    
        update_y(x, y, val, l, r, id, 0, m - 1);
    }
};

struct Vertex {
    int sum;
    Vertex *left, *right;

    Vertex(int val) : left(nullptr), right(nullptr), sum(val) {}
    Vertex(Vertex *left, Vertex *right) : left(left), right(right), sum(0) {
        if (left) sum += left->sum;
        if (right) sum += right->sum;
    }
};

struct PersistentSegmentTree {
    int n;
    vector<int> array;
    vector<Vertex*> versions;

    void build(vector<int> &array) {
        n = array.size();
        this->array = array;
        versions.push_back(build_tree(0, n - 1));
    }
    Vertex* build_tree(int l, int r) {
        if (l == r) return new Vertex(array[l]);
        return new Vertex(build_tree(l, mid), build_tree(mid + 1, r));
    }
    int answer(int L, int R, int version) {
        return answer(L, R, 0, n - 1, versions[version]);
    }
    int answer(int L, int R, int l, int r, Vertex *v) {
        if (L > R) return 0;
        if (L == l and R == r) return v->sum;
        return answer(L, min(mid, R), l, mid, v->left) + answer(max(mid + 1, L), R, mid + 1, r, v->right);    
    }
    void update(int pos, int val, int version) {
        versions.push_back(update(pos, val, 0, n - 1, versions[version]));
    } 
    Vertex* update(int pos, int val, int l, int r, Vertex *v) {
        if (l == r) return new Vertex(array[l] += val);
        if (pos <= mid) return new Vertex(update(pos, val, l, mid, v->left), v->right);
        else return new Vertex(v->left, update(pos, val, mid + 1, r, v->right));        
    }
};

struct DynamicVertex {
    int left_bound, right_bound, sum = 0;
    DynamicVertex *left_child = nullptr, *right_child = nullptr;

    DynamicVertex(int left_bound, int right_bound) {
        this->left_bound = left_bound;
        this->right_bound = right_bound;
    } 
    void extend() {
        if (!left_child and left_bound < right_bound) {
            int m = (left_bound + right_bound) >> 1;
            left_child = new DynamicVertex(left_bound, m);
            right_child = new DynamicVertex(m + 1, right_bound);
        }
    }
    void update(int pos, int val) {
        extend();
        if (left_child) {
            if (pos <= left_child->right_bound) left_child->update(pos, val);
            else right_child->update(pos, val);    
            sum = left_child->sum + right_child->sum;    
        }
        else sum += val;
    }
    int answer(int lq, int rq) {
        if (lq <= left_bound and right_bound <= rq) return sum;
        if (max(left_bound, lq) > min(right_bound, rq)) return 0;
        extend();
        return left_child->answer(lq, rq) + right_child->answer(lq, rq);    
    }
};
