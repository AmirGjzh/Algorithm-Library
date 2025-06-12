#include <bits/stdc++.h>
using namespace std;

/*--------------------------------------------------------------------------------------------------------
Breadth First Search :
it can be used for directed or undirected graphs (unweighted)
the path will be restored if it exists
Applications of BFS :

1 - Find the shortest path from a source to other vertices in an unweighted graph

2 - Find all connected components in an undirected graph in O(n + m) time: To do this, we just run BFS starting
from each vertex, except for vertices which have already been visited from previous runs. Thus, we perform normal 
BFS from each of the vertices, but do not reset the array used[] each and every time we get a new connected component, 
and the total running time will still be O(n + m)

3 - Finding a solution to a problem or a game with the least number of moves, if each state of the game can
be represented by a vertex of the graph, and the transitions from one state to the other are the edges of the graph.

4 - Finding the shortest cycle in a directed unweighted graph: Start a breadth-first search from each vertex 
As soon as we try to go from the current vertex back to the source vertex, we have found the shortest cycle containing
the source vertex. At this point we can stop the BFS, and start a new BFS from the next vertex. From all such cycles choose the shortest.

5 - Find all the edges that lie on any shortest path between a given pair of vertices (a, b)
To do this, run two breadth first searches: one from a and one from b. Let d_a[] be the array containing shortest distances 
obtained from the first BFS (from a) and d_b[] be the array containing shortest distances obtained from the second BFS from b
Now for every edge (u, v) it is easy to check whether that edge lies on any shortest path between a and b: 
the criterion is the condition d_a[u] + 1 + d_b[v] = d_a[b].

6 - Find the shortest walk of even length from a source vertex s to a target vertex t in an unweighted graph: 
For this, we must construct an auxiliary graph, whose vertices are the state (v, c), where v - the current node, 
c = 0 or c = 1 - the current parity. Any edge (u, v) of the original graph in this new column will turn into two edges 
((u, 0), (v, 1)) and ((u, 1), (v, 0)). After that we run a BFS to find the shortest walk from the starting vertex  (s, 0) 
to the end vertex (t, 0).
Note: This item uses the term "walk" rather than a "path" for a reason, as the vertices may potentially 
repeat in the found walk in order to make its length even. The problem of finding the shortest path of even length is 
NP-Complete in directed graphs, and solvable in linear time in undirected graphs, but with a much more involved approach.
Order = O(n + m)
--------------------------------------------------------------------------------------------------------*/

void BFS(int n, int s, vector<int> &d, vector<int> &p, vector<vector<int>> &g) {
    queue<int> q;
    vector<bool> used(n, false);
    q.push(s);
    p[s] = -1;
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
    for (int v = u; v != -1; v = p[u])
        path.push_back(v);
    reverse(path.begin(), path.end());
    return path;
}

/*--------------------------------------------------------------------------------------------------------
Depth First Search :
finds the lexicographical first path in the graph from a source vertex u to each vertex
Applications of DFS :

1 - Check if a vertex in a tree is an ancestor of some other vertex:
At the beginning and end of each search call we remember the entry and exit "time" of each vertex. 
Now you can find the answer for any pair of vertices (i, j) in O(1): 
vertex i is an ancestor of vertex j if and only if entry[i] < entry[j] and exit[i] > exit[j].

2 - Topological sorting:
Run a series of depth first searches so as to visit each vertex exactly once in O(n + m) time. 
The required topological ordering will be the vertices sorted in descending order of exit time.

3 - Check whether a given graph is acyclic and find cycles in a graph. (As mentioned above by 
counting back edges in every connected components).

4 - Find strongly connected components in a directed graph:
First do a topological sorting of the graph. Then transpose the graph and run another series of depth first searches in the 
order defined by the topological sort. For each DFS call the component created by it is a strongly connected component.

5 - Find bridges in an undirected graph:
First convert the given graph into a directed graph by running a series of depth first searches and 
making each edge directed as we go through it, in the direction we went. 
Second, find the strongly connected components in this directed graph. 
Bridges are the edges whose ends belong to different strongly connected components.

We perform a DFS and classify the encountered edges using the following rules:
for each edge (u, v) (u --> v):
If v is not visited:
    Tree Edge
If v is visited before u:
    Back edges - If v is an ancestor of u, then the edge
    Forward Edges - If v is a descendant of u
    Cross Edges: if v is neither an ancestor or descendant of u
Let G be an undirected graph. Then, performing a DFS upon G will classify every encountered edge as either a tree edge or back edge
Order = O(n + m)
--------------------------------------------------------------------------------------------------------*/

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
