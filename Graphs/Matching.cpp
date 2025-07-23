#include <bits/stdc++.h>
using namespace std;
using ll = long long int;

/*============================================================================================================
Bipartite Matching via Kuhn’s Algorithm

Description:
  • Finds a maximum cardinality matching in an unweighted bipartite graph
  • Uses DFS-based augmenting paths (Kuhn’s algorithm)

Applications:
  • Pairing tasks to workers, matching in bipartite networks

Notes:
  • First, check and partition the graph into two sides (0 and 1) via BFS
  • Build a reduced adjacency list from left-side vertices to right-side indices
  • Perform DFS searches to find augmenting paths, resetting the "used" marks each iteration

Order: O(V·E)
============================================================================================================*/

bool bipartite_check(const vector<vector<int>> &G, vector<int> &side) {
    int n = G.size();
    side.assign(n, -1); queue<int> q;
    for (int s = 0; s < n; s++) if (side[s] == -1) {
        q.push(s); side[s] = 0;
        while (q.size()) {
            int u = q.front(); q.pop();
            for (int v : G[u]) 
                if (side[v] == -1) {
                    side[v] = 1 - side[u];
                    q.push(v);
                } else if (side[v] == side[u]) return false;
        }
    }
    return true;
}

struct Kuhns {
    vector<bool> used;
    vector<int> matching;
    vector<vector<int>> G, adj;
    vector<pair<int, int>> match_edges;

    Kuhns(const vector<vector<int>> &adj): adj(adj) {}

    bool try_kuhn(int u) {
        if (used[u]) return false;
        used[u] = true;
        for (int v : G[u])
            if (matching[v] == -1 or try_kuhn(matching[v])) {
                matching[v] = u;
                return true;
            }
        return false;    
    }
    void max_matching() {
        int N = adj.size(), nl = 0, nr = 0;
        vector<int> side(N), lid(N, -1), rid(N, -1);
        bipartite_check(adj, side);
        for (int u = 0; u < N; u++) 
            if (side[u] == 0) lid[u] = nl++;
            else rid[u] = nr++;    
        G.assign(nl, {});
        vector<int> lo(nl), ro(nr);
        for (int u = 0; u < N; u++) {
            if (lid[u] != -1) lo[lid[u]] = u;
            if (rid[u] != -1) ro[rid[u]] = u;
        }
        for (int u = 0; u < N; u++) {
            if (side[u] == 1) continue;
            for (int v : adj[u]) G[lid[u]].push_back(rid[v]);
        }
        matching.assign(nr, -1);
        for (int u = 0; u < nl; u++) {used.assign(nl, false); try_kuhn(u);}
        for (int u = 0; u < nr; u++) if (matching[u] != -1) 
            match_edges.push_back({lo[matching[u]] + 1, ro[u] + 1});
    }
};

/*============================================================================================================
Hopcroft–Karp Maximum Bipartite Matching

Description:
  • Computes the maximum cardinality matching in an unweighted bipartite graph
  • Uses BFS to build a layered graph of shortest augmenting paths, then DFS to find multiple disjoint augmenting paths per phase

Applications:
  • Large-scale task assignment, network resource pairing, and any scenario requiring maximum bipartite pairings

Notes:
  • Pre-partition vertices into left (0…nL-1) and right (0…nR-1) sets
  • BFS sets distance labels, guiding DFS to only follow shortest paths, yielding O(√V·E) performance
  • After computing, "pairU[u]" is the matched right vertex for left "u" (or -1 if unmatched)

Order: O(√V·E)
============================================================================================================*/

struct HopcroftKarp {
    int nL, nR;
    vector<vector<int>> adj;
    vector<int> dist, pairU, pairV;
    vector<pair<int,int>> match_edges;

    HopcroftKarp(int nL, int nR)
    : nL(nL), nR(nR), adj(nL), dist(nL), pairU(nL, -1), pairV(nR, -1) {}
    void add_edge(int u, int v) {
        // u from left side (0 ... nL) and v from right side (0 ... nR)
        adj[u].push_back(v);
    }
    bool BFS() {
        queue<int> q;
        for (int u = 0; u < nL; ++u) 
            if (pairU[u] < 0) {dist[u] = 0; q.push(u);
            } else dist[u] = INT_MAX;
        bool found = false;
        while (!q.empty()) {
            int u = q.front(); q.pop();
            for (int v : adj[u]) {int pu = pairV[v];
                if (pu < 0) found = true;
                else if (dist[pu] == INT_MAX) {
                    dist[pu] = dist[u] + 1;
                    q.push(pu);
                }
            }
        }
        return found;
    }
    bool DFS(int u) {
        for (int v : adj[u]) {int pu = pairV[v];
            if (pu < 0 or (dist[pu] == dist[u] + 1 and DFS(pu))) {
                pairU[u] = v; pairV[v] = u; return true;
            }
        }
        dist[u] = INT_MAX;
        return false;
    }
    int max_matching() {
        int matching = 0;
        while (BFS()) for (int u = 0; u < nL; u++) 
            if (pairU[u] < 0 and DFS(u)) matching++;
        for (int u = 0; u < nL; u++) if (pairU[u] >= 0) 
            match_edges.emplace_back(u + 1, pairU[u] + 1);
        return matching;
    }
};

/*============================================================================================================
Hungarian Algorithm for Minimum-Cost Perfect Matching (Assignment Problem)

Description:
  • Solves the assignment problem on a bipartite graph with partitions of possibly unequal size
  • Pads the cost matrix to a square matrix of size N = max(nL, nR), filling extra cells with a large constant
  • Finds a perfect matching of minimum total cost in O(N³) time

Applications:
  • Task/worker assignments when numbers differ (unbalanced assignment)
  • Scheduling with dummy tasks or agents

Notes:
  • cost[i][j] is the cost (1-based) for matching row i to column j; pad to N×N by setting dummy cells to INF
  • Potentials u (rows) and v (columns) ensure reduced costs remain non-negative
  • Arrays p[j] and way[j] track augmenting paths during each phase
  • After computing, dummy assignments are ignored in the final output

Important:
  • Use a sufficiently large INF (e.g., 1e12) so dummy pairings are never chosen unless necessary
  • For maximum-cost assignment, negate costs before invoking solve()

Order: O(N³)
============================================================================================================*/

struct Hungarian {
    vector<int> p, way;
    vector<vector<ll>> a;  
    vector<ll> u, v, minv;
    int n, nL, nR; ll INF = 1e12 + 10;

    Hungarian(int nL, int nR, const vector<vector<ll>> &cost)
      : nL(nL), nR(nR), n(max(nL, nR)) {
        a.assign(n + 1, vector<ll>(n + 1, INF)), u.assign(n + 1, 0); 
        v.assign(n + 1, 0), minv.assign(n + 1, 0), p.assign(n + 1, 0), way.assign(n + 1, 0);
        for (int i = 1; i <= nL; ++i)
            for (int j = 1; j <= nR; ++j) a[i][j] = cost[i][j];
    }
    pair<ll, vector<int>> solve() {
        for (int i = 1; i <= n; ++i) {
            p[0] = i; int j0 = 0;
            fill(minv.begin(), minv.end(), INF);
            vector<bool> used(n + 1, false);
            do {used[j0] = true; int i0 = p[j0], j1 = 0; ll delta = INF;
                for (int j = 1; j <= n; j++) if (!used[j]) {
                    ll cur = a[i0][j] - u[i0] - v[j];
                    if (cur < minv[j]) minv[j] = cur, way[j] = j0;
                    if (minv[j] < delta) delta = minv[j], j1 = j;
                }
                for (int j = 0; j <= n; j++) 
                    if (used[j]) u[p[j]] += delta, v[j] -= delta;
                    else minv[j] -= delta;
                j0 = j1;
            } while (p[j0] != 0);
            do {int j1 = way[j0]; p[j0] = p[j1]; j0 = j1;
            } while (j0 != 0);
        }
        vector<int> assignment(nL + 1, 0); ll cost = 0;
        for (int j = 1; j <= n; j++) {int i = p[j];
            if (i >= 1 and i <= nL and j >= 1 and j <= nR) {
                assignment[i] = j;
                cost += a[i][j];
            }
        }
        return {cost, assignment};
    }
};

/*============================================================================================================
Edmonds’ Blossom Algorithm for Maximum Matching in General Graphs

Description:
  • Finds a maximum-cardinality matching in an undirected graph that may contain odd cycles
  • Contracts “blossoms” (odd-length alternating cycles) on the fly to expose augmenting paths

Applications:
  • General pairing problems beyond bipartite graphs
  • Network design, scheduling, and resource allocation with arbitrary adjacency
  • Computing minimum path covers, edge covers via matching dualities

Notes:
  • `mate[v]` holds the vertex matched with v (or –1 if unmatched)
  • BFS from each free vertex labels layers: 0 = outer, 1 = inner
  • When an edge connects two outer vertices in the same tree, we detect a blossom:
    – Compute LCA of the two endpoints in the “forest”
    – Mark all vertices on the blossom cycle and contract them by resetting their `orig[]`
  • Upon finding an unmatched vertex, we augment along the parent pointers and unwind blossoms

Complexity:
  • O(N³) in the worst-case, practical for N up to a few hundred in contests
============================================================================================================*/

struct GeneralMatching {
    int n; 
    queue<int> q;           
    vector<vector<int>> adj; 
    vector<int> mate, label, parent, orig, aux;       

    GeneralMatching(int n): n(n),
        adj(n), mate(n, -1), label(n), parent(n),
        orig(n), aux(n, -1) {
        iota(orig.begin(), orig.end(), 0);
    }
    void add_edge(int u, int v) {
        adj[u].push_back(v);
        adj[v].push_back(u);
    }
    int lca(int x, int y) {
        static int t = 0; t++;
        while (true) {
            if (x != -1) {
                if (aux[x] == t) return x;
                aux[x] = t;
                if (mate[x] == -1) x = -1;
                else x = orig[parent[mate[x]]];
            } swap(x, y);
        }
    }
    void blossom(int v, int w, int a) {
        while (orig[v] != a) {
            parent[v] = w, w = mate[v];
            if (label[w] == 1) {
                label[w] = 0; q.push(w);
            }
            orig[v] = orig[w] = a, v = parent[w];
        }
    }
    bool augment(int start) {
        fill(label.begin(), label.end(), -1);
        iota(orig.begin(), orig.end(), 0);
        while (!q.empty()) q.pop();
        label[start] = 0; parent[start] = -1; q.push(start);
        while (!q.empty()) {
            int v = q.front(); q.pop();
            for (int x : adj[v])
                if (label[x] == -1) {
                    label[x] = 1, parent[x] = v;
                    if (mate[x] == -1) {int j = x;
                        while (j != -1) {
                            int p = parent[j], nxt = (p == -1 ? -1 : mate[p]);
                            mate[j] = p; mate[p] = j; j = nxt;
                        }
                        return true;
                    }
                    label[mate[x]] = 0; q.push(mate[x]);
                } else if (label[x] == 0 and orig[v] != orig[x] and mate[v] != x) {
                    int a = lca(orig[v], orig[x]);
                    blossom(x, v, a); blossom(v, x, a);
                }
        }
        return false;
    }
    int max_matching() {
        int match_count = 0;
        for (int v = 0; v < n; v++) 
            if (mate[v] == -1) for (int u : adj[v]) if (mate[u] == -1) {
                mate[v] = u; mate[u] = v; match_count++;
                break;
            }
        for (int v = 0; v < n; v++) 
            if (mate[v] == -1 and augment(v)) match_count++;
        return match_count;
    }
};

/*============================================================================================================
Kőnig’s Theorem: Min Vertex Cover = Max Matching in Bipartite Graphs

Statement:
  • In any bipartite graph, the size of a maximum matching equals the size of a minimum vertex cover

Definitions:
  – Matching: A set of edges without common vertices
  – Vertex cover: A set of vertices such that every edge has at least one endpoint in the set

Implications:
  • After computing a maximum matching (e.g., via Kuhn’s or Hopcroft–Karp), one can derive a minimum vertex cover:
    1. Find all free (unmatched) vertices on the left side
    2. Perform alternating BFS/DFS over “non‑matching” edges from left-to-right and “matching” edges from right-to-left
    3. Let L+ be the reachable left vertices, and R+ be the reachable right vertices
    4. Then the minimum vertex cover is:
      (All left vertices \ L+)  ∪  (R+)

Applications:
  • Bipartite scheduling: minimize resources covering all pairings
  • Network security: minimal guard placement

Use this theorem to extract both the matching and the corresponding minimum vertex cover in O(V+E) after your matching algorithm
============================================================================================================*/
