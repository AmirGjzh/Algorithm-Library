#include <bits/stdc++.h>
using namespace std;

/*============================================================================================================
Breadth‑First Search (BFS)

Description:
  • Works on both directed and undirected unweighted graphs
  • Can reconstruct the shortest path from a source to any reachable vertex

Applications:
  1. Shortest paths in an unweighted graph:
    – Run once from source s to compute dist[] and parent[] in O(n + m)

  2. Connected components in an undirected graph:
    – Reuse the same visited[] array across runs
    – Start BFS from each unvisited vertex, total time still O(n + m)

  3. Minimum‑move solutions in state‑space search:
    – Model each state as a vertex, transitions as edges, then BFS for fewest moves

  4. Shortest directed cycle through a given vertex:
    – For each vertex u:
      • Run BFS from u
      • As soon as you see an edge back to u, you’ve found the shortest cycle containing u
  
  5. Identifying edges on any shortest path between a and b:
    – Run BFS from a → d_a[]
    – Run BFS from b → d_b[]
    – Edge (u, v) lies on some shortest a → b iff  
      • d_a[u] + 1 + d_b[v] == d_a[b]

  6. Shortest even‑length walk (undirected):
    – Build an “auxiliary” graph on states (v, parity)
    – Replace each original edge (u, v) with:
      • ((u, 0) → (v, 1)) and ((u, 1) → (v, 0))
    – BFS from (s, 0) to (t, 0) finds shortest even‑length walk
      • (Vertices may repeat; finding a simple even‑length path is NP‑complete.)
      
Order: O(n + m)
============================================================================================================*/

void BFS(int n, int s, vector<int> &d, vector<int> &p, vector<vector<int>> &g) {
    queue<int> q;
    d.assign(n, 1e9 + 10);
    p.assign(n, -1);
    vector<bool> used(n, false);
    q.push(s);
    used[s] = true;
    while (q.size()) {
        int u = q.front();
        q.pop();
        for (int v : g[u]) 
            if (!used[v]) {
                used[v] = true;
                q.push(v);
                d[v] = d[u] + 1;
                p[v] = u;
            }
    }
}

vector<int> BFS_path(int u, vector<int> &p) {
    vector<int> path;
    for (int v = u; v != -1; v = p[v])
        path.push_back(v);
    reverse(path.begin(), path.end());
    return path;
}

/*============================================================================================================
Depth‑First Search (DFS)

Description:
  • Explores as deeply as possible before backtracking
  • Records entry and exit “times” for each vertex
  
Applications:
  1. Ancestor queries in a tree:
    – Track entry[u], exit[u] during DFS
    – u is ancestor of v iff entry[u] < entry[v] and exit[u] > exit[v]

  2. Topological sorting (directed acyclic graph):
    – After a full DFS, sort vertices in descending order of exit time

  3. Cycle detection in directed graphs:
    – Back edges (to a gray ancestor) indicate a cycle

  4. Strongly connected components (Kosaraju’s algorithm):
    – Run DFS to compute exit times
    – Transpose graph
    – Process vertices in order of decreasing exit time, each DFS marks one SCC

Edge classification (in undirected graph G):
  – Tree edges:    to an unvisited vertex
  – Back edges:    to a gray ancestor
  – Forward edges: to a black descendant (in directed graphs)
  – Cross edges:   to a black non-descendant (in directed graphs)

Order: O(n + m)
============================================================================================================*/

int n, timer = 0;
vector<vector<int>> g;
vector<int> in(n), out(n), color(n, 0);

void DFS(int u) {
    in[u] = timer++;
    color[u] = 1;
    for (int v : g[u]) 
        if (color[v] == 0) 
            DFS(v);
    color[u] = 2;
    out[u] = timer++;    
}
