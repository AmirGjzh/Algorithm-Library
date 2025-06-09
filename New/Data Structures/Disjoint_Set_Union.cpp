#include <bits/stdc++.h>
using namespace std;
const int N = 1e6 + 10;

/*--------------------------------------------------------------------------------------------------------
DSU :
path compression => makes the order better (Log(n) and 1 with size) but we loose parents of some datas
using size => makes the order better (Log(n) and 1 with path compression)
not using size, will give use better flexibility (we can union in a way we want)
DSU allows you to easily store additional information in the sets
--------------------------------------------------------------------------------------------------------*/

struct Node {
    int parent;
    int size;
};

struct DisjointSetUnion {
    vector<Node> node;

    void build(int n) {
        node.resize(n);
        for (int i = 0; i < n; i++)
            make_set(i);
    }

    void make_set(int a) {
        node[a].parent = a;
        node[a].size = 1;
    }

    int find_set(int a) {
        if (a != node[a].parent) 
            node[a].parent = find_set(node[a].parent);
        return node[a].parent;
    }

    void union_set(int a, int b) {
        a = find_set(a);
        b = find_set(b);
        if (a != b) 
            node[b].parent = a;
    }

    void union_set_balance(int a, int b) {
        a = find_set(a);
        b = find_set(b);
        if (a != b) {
            if (node[a].size < node[b].size)
                swap(a, b);
            node[b].parent = a;
            node[a].size += node[b].size;    
        }
    }
};
