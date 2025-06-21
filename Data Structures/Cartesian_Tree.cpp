#include <bits/stdc++.h>
using namespace std;

/*============================================================================================================
Implicit Treap (Cartesian Tree)

Description:
  • A randomized BST + heap ("treap") that supports array-like operations:
    – point insert/delete at any position
    – range queries (sum, min, max, etc.)
    – range updates (add, set, reverse, etc.)
  • Powered by split/merge + lazy propagation for efficiency

Core Concepts:
  – Keys are implicit: defined by in-order index (0-based), not stored explicitly
  – Random priority ensures expected O(log(n)) time per operation via heap property
  – Maintains subtree metadata (size, sum) and lazy flags (assign, reverse)

Operations:
  * split(root, left, right, key):
    – Divides treap into [0..key-1] and [key..n-1], updating size and pushing laziness
  * merge(root, left, right):
    – Combines two treaps into one; left < right by indices, pushes laziness before merging
  * insert(pos, value): split at pos → merge(left, new node) → merge with right
  * remove(pos): split at pos; split right by 1 → discard middle → merge rest
  * range update (add x to [L..R]): split into a, b, c → mark b's lazy flag → merge back
  * range query: split into a, b, c → answer b.sum → merge back

Lazy Handling (`push`):
  • Reverse flag: swaps children and propagates toggles downstream
  • Add-update flag: accumulates in `lazy`, marks node, later applied in `pull`

Complexity:
  • All operations (split, merge, insert, delete, update, query) run in expected O(log(n)) time

Use Cases:
  • Dynamic sequence manipulation (e.g., ropes, substring reversals, interval updates/queries)

Notes:
  • Works on implicit indices, keys = size(left subtree).
  • Always `push` before traversing `left`/`right` to ensure correct metadata
  • `pull` updates `size` and `sum` based on children after modifications
============================================================================================================*/

struct Node {
    int val;
    int sum;
    int lazy; 
    bool mark, reverse;

    int size, priority;
    Node *left, *right, *par;
    Node(int val) {
        mark = reverse = false;
        this->val = sum = val;
        lazy = 0, size = 1, priority = rand();
        left = right = par = nullptr;
    }
};

struct Treap {
    Node *root;

    void build(vector<int> &a) {
        root = nullptr;
        for (int i : a) 
            merge(root, root, new Node(i));
    }

    int size(Node *root) {
        return root ? root->size : 0;
    }

    int sum(Node *root) {
        return root ? root->sum : 0;
    }

    void push(Node *root) {
        if (!root) return;
        if (root->reverse) {
            root->reverse = false;
            swap(root->left, root->right);
            if (root->left) root->left->reverse ^= true;
            if (root->right) root->right->reverse ^= true;
        }
        if (root->left and root->mark) {
            root->left->lazy += root->lazy;
            root->left->mark = true;
        }
        if (root->right and root->mark) {
            root->right->lazy += root->lazy;
            root->right->mark = true;
        }
        if (root->mark) {
            root->sum += size(root) * root->lazy;
            root->lazy = 0, root->mark = false;
        }
    }
    
    void pull(Node *root) {
        if (!root) return;
        push(root->left), push(root->right);
        root->size = size(root->left) + size(root->right) + 1;
        root->sum = sum(root->left) + sum(root->right) + root->val;
    }

    void split(Node *root, Node *&left, Node *&right, int key) {
        if (!root)
            return void(left = right = nullptr);
        push(root);    
        if (size(root->left) < key) {
            if (root->right) root->right->par = nullptr;
            split(root->right, root->right, right, key - size(root->left) - 1);
            if (root->right) root->right->par = root;
            left = root;
        }
        else {
            if (root->left) root->left->par = nullptr;
            split(root->left, left, root->left, key);
            if (root->left) root->left->par = root;
            right = root;
        }
        pull(root);    
    }

    void merge(Node *&root, Node *left, Node *right) {
        push(left), push(right);
        if (!left or !right)
            return void(root = left ? left : right);
        else if (left->priority < right->priority) {
            if (left->right) left->right->par = nullptr;
            merge(left->right, left->right, right);
            if (left->right) left->right->par = left;
            root = left;
        }
        else {
            if (right->left) right->left->par = nullptr;
            merge(right->left, left, right->left);
            if (right->left) right->left->par = right;
            root = right; 
        }
        pull(root);           
    }

    void insert(int ind, int val) {
        Node *l, *r;
        split(root, l, r, ind);   
        Node *node = new Node(val);
        merge(l, l, node);
        merge(root, l, r);
    }

    void remove(int ind) {
        Node *l, *r;
        split(root, l, r, ind);
        split(r, root, r, 1);
        merge(root, l, r);
    }

    void update(int L, int R, int val) {
        Node *a, *b, *c;
        split(root, a, b, L);
        split(b, b, c, R - L + 1);
        b->lazy += val, b->mark = true;
        merge(root, a, b);
        merge(root, root, c);
    }

    int answer(int L, int R) {
        Node *a, *b, *c;
        split(root, a, b, L);
        split(b, b, c, R - L + 1);
        int res = b->sum;
        merge(root, a, b);
        merge(root, root, c);
        return res;
    }

    void output(Node *root) {
        if (!root) return;
        push(root);
        output(root->left);
        cout << root->val;
        output(root->right);
    }
};
