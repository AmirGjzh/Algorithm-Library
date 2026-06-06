#include <bits/stdc++.h>
using namespace std;

/*============================================================================================================
Disjoint Set Union (DSU / Union-Find)

Description:
  • Maintains a collection of disjoint (non-overlapping) sets over elements [0..n-1]
  • Each element belongs to exactly one set, represented by a single root (leader)
  • Supports efficient merging of sets and membership queries

Core Operations:
  1. make_set(x)
    • Initializes element x as a standalone set
    • Parent of x is itself, size initialized to 1

  2. find_set(x)
    • Returns the representative (root) of the set containing x
    • Uses path compression to flatten the tree for faster future queries

  3. union_set(a, b)
    • Merges the sets containing a and b
    • Attaches the root of b under the root of a
    • Does not use size or rank heuristics

  4. union_set_balance(a, b)
    • Merges sets using union by size
    • Smaller set is attached under the larger one
    • Maintains size information for balanced trees

Data Representation:
  • Each element stores:
    – parent: pointer to its parent (or itself if root)
    – size: number of elements in the set (valid when using balanced unions)

Optimizations:
  • Path Compression:
    – Applied in find_set
    – Ensures near-constant amortized time per operation

Time Complexity:
  • make_set: O(1)
  • find_set: Amortized O(α(n)) with path compression
  • union_set / union_set_balance: Amortized O(α(n))

Applications:
  • Connected component detection in graphs
  • Kruskal’s Minimum Spanning Tree algorithm
  • Cycle detection in undirected graphs
  • Network connectivity and clustering problems

Notes:
  • Elements are indexed from 0 to n-1
  • Representative elements are not fixed and may change after unions
============================================================================================================*/

struct Node {
    int parent, size;
};

struct DisjointSetUnion {
    vector<Node> node;

    explicit DisjointSetUnion(const int n): node(n) {
        for (int i = 0; i < n; i++) make_set(i);
    }
    void make_set(const int a) {
        node[a].parent = a, node[a].size = 1;
    }
    int find_set(const int a) {
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
