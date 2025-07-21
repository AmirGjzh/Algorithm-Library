#include <bits/stdc++.h>
using namespace std;

/*============================================================================================================
Disjoint Set Union (DSU)

Description:
  DSU is a data structure that keeps track of elements partitioned into disjoint sets
  It supports two main operations efficiently:
    - find_set(a): Find the representative (root) of the set containing 'a' with path compression
    - union_set(a, b): Merge the sets containing 'a' and 'b'

Features and Variants:
  • Path Compression: Flattens the tree to speed up future queries, amortized nearly O(1)
  • Union by Size (or Rank): Always attach smaller tree under the root of the larger tree to keep trees shallow
  • Without union by size: allows flexible union strategies but can degrade performance

Structure:
  - Node struct stores 'parent' and 'size' (size of the set for union balancing)
  - DisjointSetUnion class encapsulates DSU functionality with convenient methods

Applications:
  • Connected components in graphs
  • Kruskal's MST algorithm
  • Detecting cycles in undirected graphs
  • Network connectivity
============================================================================================================*/

struct Node {
    int parent, size;
};

struct DisjointSetUnion {
    vector<Node> node;

    DisjointSetUnion(const int n): node(n) {
        for (int i = 0; i < n; i++) make_set(i);
    }
    void make_set(int a) {
        node[a].parent = a, node[a].size = 1;
    }
    int find_set(int a) {
        if (a != node[a].parent) node[a].parent = find_set(node[a].parent);
        return node[a].parent;
    }
    void union_set(int a, int b) {
        a = find_set(a), b = find_set(b);
        if (a != b) node[b].parent = a;
    }
    void union_set_balance(int a, int b) {
        a = find_set(a), b = find_set(b);
        if (a != b) {
            if (node[a].size < node[b].size) swap(a, b);
            node[b].parent = a;
            node[a].size += node[b].size;    
        }
    }
};
