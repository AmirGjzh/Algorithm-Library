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
    d.assign(n, 1e9 + 10), p.assign(n, -1);
    vector<bool> used(n, false);
    q.push(s); used[s] = true;
    while (q.size()) {
        int u = q.front();
        q.pop();
        for (int v : g[u]) 
            if (!used[v]) {
				q.push(v);
                used[v] = true, d[v] = d[u] + 1, p[v] = u;
            }
    }
}

vector<int> BFS_path(int u, vector<int> &p) {
    vector<int> path;
    for (int v = u; v != -1; v = p[v]) path.push_back(v);
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
    for (int v : g[u]) if (color[v] == 0) DFS(v);
    color[u] = 2;
    out[u] = timer++;    
}

/*============================================================================================================
Cycle Detection in Graphs (Directed & Undirected)

Description:
  • Provides methods to detect a cycle in both directed and undirected graphs
  • Uses depth-first search with vertex coloring for directed graphs
  • Uses parent-tracking DFS for undirected graphs to avoid trivial backtracking

Applications:
  • Detecting infinite loops or deadlocks in dependency graphs
  • Verifying DAG property in scheduling or prerequisite structures
  • Checking connectivity constraints where cycles are forbidden (e.g., tree validation)

Notes:
  • `col[u]` coloring: 0 = unvisited, 1 = in-stack (active), 2 = fully explored
  • In directed DFS, encountering a neighbor with `col[v] == 1` signals a back-edge (cycle)
  • In undirected DFS, skip the edge back to parent, encountering any other visited neighbor signals a cycle
  • On detection, records `cycle_start` and `cycle_end`, then reconstructs the cycle path via the `par` array

Order:
  • Directed DFS: O(n + m)
  • Undirected DFS: O(n + m)
============================================================================================================*/

struct FindCycle {
	vector<int> col, par;
	vector<vector<int>> g;
	int n, cycle_start, cycle_end;

	bool DFS_directed(int u) {
		col[u] = 1;
		for (int v : g[u]) 
			if (col[v] == 0) {par[v] = u; if (DFS_directed(v)) return true;}
			else if (col[v] == 1) {cycle_end = u, cycle_start = v; return true;}
		col[u] = 2;
		return false;	
	}
	bool DFS_undirected(int u, int p) {
		col[u] = 1;
		for (int v : g[u]) {
			if (v == p) continue;
			if (col[v] == 1) {cycle_end = u, cycle_start = v;return true;}
			par[v] = u;
			if (DFS_undirected(v, u)) return true;
		}
		return false;
	}
	void find_cycle() {
		col.assign(n, 0), par.assign(n, -1);
		cycle_start = cycle_end = -1;
		for (int u = 0; u < n; u++) if (col[u] == 0 and DFS_directed(u)) break;
		if (cycle_start == -1) cout << "Acyclic";
		else {
			vector<int> cycle;
			cycle.push_back(cycle_start);
			for (int u = cycle_end; u != cycle_start; u = par[u]) cycle.push_back(u);
			reverse(cycle.begin(), cycle.end());
		}	
	}
};

/*============================================================================================================
Eulerian Path & Cycle (Recursive Hierholzer’s Algorithm)

Description:
  • Recursively constructs an Eulerian path or cycle in an undirected multigraph (supports multiple edges and loops)
  • Visits each edge exactly once, appending vertices post‐order to build the tour in reverse

Applications:
  • Route inspection tasks (e.g., garbage collection, street sweeping)
  • Puzzle solving (drawing figures with a single stroke)
  • Network traversal where each link must be traversed exactly one time

Notes:
  • Graph must be connected when ignoring isolated vertices
  • Exactly 0 or 2 vertices of odd degree → cycle or open trail respectively
  • Start must be one of the odd nodes, for tour, or any node for cycle
  • Each undirected edge is assigned one unique ID, marking it prevents revisiting
  • Resulting `tour` vector is generated in reverse—remember to reverse before use
  • Recursive depth = length of tour, for very large graphs consider the iterative version

Order: O(n + m)
============================================================================================================*/

struct EulerTour {
    int n, m;                               
    vector<bool> used;                          
    vector<int> ptr, tour;                        
    vector<vector<pair<int,int>>> adj;    

    void add_edge(int u, int v, int eid) {
        adj[u].emplace_back(v, eid);
        adj[v].emplace_back(u, eid);
    }
    void DFS(int u) {
        while (ptr[u] < (int) adj[u].size()) {
            auto [v, eid] = adj[u][ptr[u]++];
            if (!used[eid]) {used[eid] = true; DFS(v);}
        }
        tour.push_back(u);
    }
    vector<int> get_tour(int start) {
        ptr.assign(n, 0), used.assign(m, false);
        DFS(start);
        reverse(tour.begin(), tour.end());
        return tour;
    }
};