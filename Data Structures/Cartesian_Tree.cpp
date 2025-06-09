#include <bits/stdc++.h>
using namespace std;

/*--------------------------------------------------------------------------------------------------------
Treap :
A randomized binary search tree
It can do all the things segment tree does
like point and range queries
plus some more beautiful queries like inverting a segment
It can cut and paste a segment
It can insert or remove an elemnt
All we need :)
we do all our works with split and merge function
pull function is the same as combine in segtree
push function is the same as push in lazysegtree
but we have to be carefull with pushing more than 1 lazy values (order of them is important)
Orders :
split => O(Log(n))
merge => O(Log(n))
--------------------------------------------------------------------------------------------------------*/

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
