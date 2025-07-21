#include <bits/stdc++.h>
using namespace std;
using ll = long long int;

/*============================================================================================================
Edmonds–Karp Maximum Flow (Ford–Fulkerson with BFS)

Description:
  • Computes the maximum flow from source s to sink t in a directed graph
  • Uses repeated BFS to find the shortest (fewest‑edge) augmenting path in the residual network
  • Each BFS identifies a bottleneck capacity; we augment along that path and update residual capacities

Applications:
  • Network routing & bandwidth allocation
  • Bipartite matching (via reduction to max‑flow)
  • Circulation with lower/upper bounds and project‐selection problems

Complexity:
  • O(V·E²) worst‐case:
    – Each BFS is O(E) to find an augmenting path
    – There can be up to O(V·E) augmentations in the Ford–Fulkerson method

Implementation Details:
  • G[u] holds neighbors for both forward and backward edges
  • cap[u][v]: current residual capacity from u→v
  • flow[u][v]: total flow sent along u→v (negative on reverse edge)

Edge‐Case Handling:
  • Multiple parallel edges (u→v added twice):  
    – cap[u][v] is incremented, not overwritten, the adjacency list is added only on the first call, avoiding duplicates
  • Opposite‐direction edges (u→v and v→u):  
    – Both appear in G[u] and G[v], They each maintain independent capacities, so you can send flow both ways correctly
  • Self‐loops (u==v):  
    – We insert u into G[u] once, BFS will ignore it (par[u] already set), so it has no effect on augmentations
============================================================================================================*/

struct EdmondsKarp {
    int n;
    vector<vector<int>> G;
    vector<vector<ll>> cap, flow;

    EdmondsKarp(int n): n(n) {
        G.assign(n, {});
        cap.assign(n, vector<ll>(n, 0));
        flow.assign(n, vector<ll>(n, 0));
    }
    void add_edge(int u, int v, ll c) {
        if (cap[u][v] == 0 and cap[v][u] == 0) 
            G[u].push_back(v), G[v].push_back(u);
        cap[u][v] += c;
    }
    ll max_flow(int S, int T) {
        ll total = 0; vector<int> par(n);
        while (true) {
            fill(par.begin(), par.end(), -1);
            par[S] = -2; queue<pair<ll, ll>> q;
            q.push({S, LLONG_MAX}); ll new_flow = 0;
            while (q.size()) {
                auto [u, f] = q.front(); q.pop();
                for (int v : G[u]) 
                    if (par[v] < 0 and cap[u][v] > 0) {
                        par[v] = u; ll mf = min(f, cap[u][v]);
                        if (v == T) {new_flow = mf; break;}
                        q.push({v, mf});
                    }
                if (new_flow) break;    
            }
            if (!new_flow) break;
            total += new_flow;
            int cur = T;
            while (cur != S) {
                cap[par[cur]][cur] -= new_flow;
                cap[cur][par[cur]] += new_flow;
                flow[par[cur]][cur] += new_flow;
                flow[cur][par[cur]] -= new_flow;
                cur = par[cur];
            }
        }
        return total;
    }
};

/*============================================================================================================
Push–Relabel Maximum Flow (FIFO Preflow–Push)

Description:
  • Computes the maximum flow from source S to sink T in a directed graph using a preflow‑push approach

Key Data Structures:
  – cap[u][v]: original capacity of edge u→v  
  – flow[u][v]: current net flow on u→v (can be negative on the reverse direction)  
  – excess[u]: net inflow minus outflow at u  
  – height[u]: distance label (an integer “height”)  
  – seen[u]: index of next neighbor to try during discharge  
  – active: FIFO queue of vertices (other than S,T) with excess > 0  

Complexity:
  • O(V²·E) worst‑case  
    – Each discharge scans up to E edges; each relabel raises height[u] by at least 1  
    – Height can increase at most 2V–1 times per vertex, so total relabels ≤ O(V²)
    – Pushes and relabels combined give O(V²·E)

Behavior on Special Cases:
  • Parallel edges (u→v added multiple times):  
    Each call to add_edge(u,v,c) simply increments cap[u][v] by c
    During pushes, residual capacity is cap[u][v]−flow[u][v], so the algorithm naturally treats  
    the combined parallel edges as one aggregated capacity
  • Opposite‐direction edges (u→v and v→u both present):  
    Those are stored in cap[u][v] and cap[v][u] separately, A push from u→v increases flow[u][v]  
    and decreases flow[v][u], and vice versa, correctly handling a true bidirectional network  
  • Self‐loops (u→u): 
    Self‑loops contribute cap[u][u], but since height[u] = height[u] + 1 is never true,  
    no push will ever occur on a self‑loop—so they safely remain unused and do not corrupt labeling
============================================================================================================*/

struct PushRelabel {
    int n;
    vector<ll> excess;
    queue<int> active;
    vector<int> height, seen;
    vector<vector<ll>> cap, flow;

    PushRelabel(int n): n(n) {
        excess.assign(n, 0);
        cap.assign(n, vector<ll>(n, 0));
        flow.assign(n, vector<ll>(n, 0));
        seen.assign(n, 0), height.assign(n, 0);
    }
    void add_edge(int u, int v, ll c) {
        cap[u][v] += c;
    }
    void push(int u, int v, int T) {
        ll send = min(excess[u], cap[u][v] - flow[u][v]);
        if (send <= 0 or height[u] != height[v] + 1) return;
        excess[u] -= send, excess[v] += send;
        flow[u][v] += send, flow[v][u] -= send;
        if (excess[v] == send and v != T) active.push(v);
    }
    void relabel(int u) {
        ll mh = LLONG_MAX;
        for (int v = 0; v < n; v++) 
            if (cap[u][v] - flow[u][v] > 0) 
                mh = min(mh, (ll) height[v]);
        if (mh < LLONG_MAX) height[u] = (int) mh + 1;        
    }
    void discharge(int u, int T) {
        while (excess[u] > 0) 
            if (seen[u] < n) {int v = seen[u]++; push(u, v, T);}
            else {relabel(u); seen[u] = 0;}
    }
    ll max_flow(int S, int T) {
        height[S] = n;
        for (int u = 0; u < n; u++) 
            if (cap[S][u] > 0) {
                flow[S][u] = cap[S][u], flow[u][S] = -flow[S][u];
                excess[u] = flow[S][u], excess[S] -= flow[S][u];
                if (u != T) active.push(u);
            }
        while (active.size()) {
            int u = active.front(); active.pop();
            if (u != S and u != T) discharge(u, T);
        }
        return excess[T];
    }
};


/*============================================================================================================
Dinic's Algorithm for Maximum Flow

Description:
  • Computes the maximum flow from source S to sink T using level graphs and blocking flows

Expected Performance:
  • In practice, often O(E·√V) on sparse or random graphs
  • For unit capacities or bipartite matching, guarantees O(√V·E)
  • Worst‐case O(V²·E), but pathological examples are uncommon in contest use

Handling Special Cases:
  1. Parallel edges (multiple u→v calls):
    Each add_edge(u,v,cap) pushes a new Edge struct onto G[u] and a reverse zero‐cap
    entry on G[v], The algorithm naturally aggregates these parallel edges since level‐
    and capacity‐based decisions treat each entry independently, effectively summing
    capacities
  2. Bidirectional edges (both u→v and v→u):**
    Internally these become two distinct forward edges plus their zero‐cap reverse half:
      - The `add_edge(u,v,cap)` call creates cap on u→v and a zero‐cap reverse v→
      - A separate `add_edge(v,u,cap2)` call likewise sets up cap2 on v→u and zero‐cap u→v
    Residual capacities on each direction evolve correctly, allowing flow to route
    back and forth
  3. Self‑loops (u→u):
    A self‐loop entry is added to G[u] with nonzero capm However, since Dinic’s BFS
    advances only to neighbors with level +1, and level[u] never equals level[u]+1,
    no DFS can ever push flow along a self‐loop, It remains inert and harmless
============================================================================================================*/

struct Dinic {
    struct Edge {
        int to, rev; ll cap;
    };

    int n, S, T;
    vector<int> level, ptr;
    vector<vector<Edge>> G;
    vector<vector<ll>> flow;
    vector<tuple<int, int, ll>> edges;

    Dinic(int n, int S, int T): n(n), S(S), T(T) {
        flow.assign(n, vector<ll>(n, 0));
        G.assign(n, {}); level.assign(n, 0); ptr.assign(n, 0);
    }
    void add_edge(int u, int v, ll cap) {
        G[u].push_back({v, (int) G[v].size(), cap});
        G[v].push_back({u, (int) G[u].size() - 1, 0});
        edges.emplace_back(u, v, cap);
    }
    bool BFS() {
        fill(level.begin(), level.end(), -1);
        level[S] = 0; queue<int> q; q.push(S);
        while (q.size()) {
            int u = q.front(); q.pop();
            for (auto &e : G[u]) 
                if (level[e.to] < 0 and e.cap > 0) {
                    level[e.to] = level[u] + 1;
                    q.push(e.to);
                }
        }
        return level[T] >= 0;
    }
    ll DFS(int u, ll pushed) {
        if (!pushed or u == T) return pushed;
        for (int &cid = ptr[u]; cid < G[u].size(); cid++) {
            auto &e = G[u][cid];
            if (level[e.to] != level[u] + 1 or e.cap == 0) continue;
            ll tr = DFS(e.to, min(pushed, e.cap));
            if (tr) {
                e.cap -= tr; G[e.to][e.rev].cap += tr;
                flow[u][e.to] += tr; flow[e.to][u] -= tr;
                return tr;
            }
        }
        return 0;
    } 
    ll max_flow() {
        ll flow = 0;
        while (BFS()) {
            fill(ptr.begin(), ptr.end(), 0);
            while (ll pushed = DFS(S, LLONG_MAX)) flow += pushed;
        }
        return flow;
    }
};

/*============================================================================================================
Edmonds–Karp with Lower Bounds & Flow Recovery (FeasibleFlow)

Description:
  • Computes a feasible circulation on a directed graph with per‐edge lower (L) and upper (C) bounds

Applications:
  • Circulation with demands (supply/demand networks)
  • Scheduling/resource allocation with per‐task minimum and maximum process rates
  • Any scenario requiring both lower‐ and upper‐bounded flows

Complexity:
  • Standard Edmonds–Karp on a graph with N = original nodes + 2 auxiliary (SS, TT):
    – O(N·E²) worst‐case
    – Typically performs well for modest N, E (≤ several thousands)

Handling Special Cases:
  1. Multiple parallel edges u→v
    • Successive `add_edge(u,v,c)` calls accumulate `cap[u][v] += c`
    • In the adjacency `G[u]` only a single entry is created on the first such call, 
       but total capacity sums correctly
  2. Bidirectional edges u→v and v→u
    • Both directions are stored independently in G and cap[][], allowing true two‐
      way flows, A call to add_edge(u,v,c1) and add_edge(v,u,c2) yields cap[u][v]=c1, cap[v][u]=c2
  3. Self‑loops (u→u)  
    • Added to G[u] once, however, BFS only advances along edges with positive residual 
      *and* parent[v] unset, Since for a self‐loop v==u and parent[u] is already set, 
      no augmenting path uses it, It remains inert—harmless but counted in adjacency
============================================================================================================*/

struct FeasibleFlow {
    int n;
    vector<vector<int>> G;
    vector<vector<ll>> cap, flow;

    FeasibleFlow(int n): n(n) {
        G.assign(n, {});
        cap.assign(n, vector<ll>(n, 0));
        flow.assign(n, vector<ll>(n, 0));
    }
    void add_edge(int u, int v, ll c) {
        if (cap[u][v] == 0 and cap[v][u] == 0) {
            G[u].push_back(v);
            G[v].push_back(u);
        }
        cap[u][v] += c;
    }
    ll BFS(int S, int T, vector<int> &par) {
        fill(par.begin(), par.end(), -1);
        par[S] = -2; queue<pair<int, ll>> q;
        q.push({S, LLONG_MAX});
        while (q.size()) {
            auto [u, f] = q.front(); q.pop();
            for (int v : G[u]) 
                if (par[v] == -1 and cap[u][v] > 0) {
                    par[v] = u; ll newf = min(f, cap[u][v]);
                    if (v == T) return newf;
                    q.push({v, newf});
                }
        }
        return 0;
    }
    ll max_flow(int S, int T) {
        vector<int> par(n);
        ll flow_sum = 0, newf;
        while ((newf = BFS(S, T, par)) > 0) {
            flow_sum += newf;
            int cur = T; while (cur != S) {
                cap[par[cur]][cur] -= newf;
                cap[cur][par[cur]] += newf;
                flow[par[cur]][cur] += newf;
                flow[cur][par[cur]] -= newf;
                cur = par[cur];
            }
        }
        return flow_sum;
    }
    ll solve(int m, int S, int T) {
        vector<ll> balance(n, 0);
        vector<tuple<int,int,ll,ll>> edges;
        for (int i = 0; i < m; i++) {
            int u, v; ll L, C;
            cin >> u >> v >> L >> C;
            edges.emplace_back(--u, --v, L, C);
            balance[u] -= L; balance[v] += L;
        }
        int SS = n, TT = n + 1; ll totalSupply = 0;
        FeasibleFlow ff(n + 2);
        for (int i = 0; i < n; i++) 
            if (balance[i] > 0) {
                ff.add_edge(SS, i, balance[i]);
                totalSupply += balance[i];
            } else if (balance[i] < 0) ff.add_edge(i, TT, -balance[i]);
        for (auto &e: edges) {
            int u, v; ll L, C; tie(u, v, L, C) = e;
            ff.add_edge(u, v, C - L);
        }
        ff.add_edge(T, S, LLONG_MAX);
        ll fl = ff.max_flow(SS, TT);
        if (fl != totalSupply) {
            cout << "No feasible flow exists\n";
            return - 1;
        }
        return ff.max_flow(S, T);
    }
};

/*============================================================================================================
Min–Cost Max–Flow with Lower Bounds (Successive Shortest Path + Johnson Potentials)

Problem Statement:
  • Given a directed graph with n nodes (0…n−1), you must send up to K units of flow
    from a source S to a sink T, minimizing total cost
  • Each original edge (u→v) has:
    – a mandatory minimum flow `low` (you must send at least this many units),
    – a capacity `up` (cannot exceed this),
    – and a per-unit cost `cost`

Technique Overview: 
  2) Feasibility Check
    – Introduce auxiliary super‑source SS = n and super‑sink TT = n+1
    – For every i with demand[i] > 0, add SS→i of cap = demand[i], cost = 0  
    – For every i with demand[i] < 0, add i→TT of cap = −demand[i], cost = 0  
    – Run a basic MinCostFlow(SS, TT) for total demand = sum of all positive demands 
      If max‐flow < total demand, no circulatio* can satisfy all lower bounds  
  3) Maximization Phase
    – If feasible, add an infinite‐capacity, zero‐cost back‑edge T→S to allow further flow  
    – Run MinCostFlow(S, T, K) to push up to K additional units  
    – The combined flows respect all original `[low, up]` constraints, minimize cost,  
      and never violate capacity or lower‑bound requirements

Time & Memory:
  • Feasibility check: one MCF run on n+2 nodes, E+|demands| edges → O(E.log(V))
  • Main solve: another MCF run up to K units → O(min(K, totalAugments)·E·log(V))  
  • Memory: O((n+2) + E) for adjacency and reverse‑edge lists, plus O(n) for potentials, distances

Edge‐Case Handling:
  1. Multiple parallel edges (u→v)
    – Each add_edge(u,v,low,up,cost) call appends a fresh entry to `G[u]`, `G[v]`
    – Flows on them are tracked independently, their capacities simply coexist
  2. Opposite‐direction edges (v→u and u→v)
    – As above, each direction is stored distinctly, Residual capacity and cost correctly 
       handle pushing flow in either direction without interference
  3. Self‑loops (u→u)
    – They appear in `G[u]` and carry zero reduced‑cost cycles in practice
    – Dijkstra/potentials skip them because residual capacity > 0 but `u==v` → no distance gain
    – They do not affect shortest‐path distances or net flow
============================================================================================================*/

struct MinCostFlow {
    struct Edge {
        int to, rev; ll cap, cost;
    };

    int n;
    vector<ll> dis, pot;
    vector<vector<Edge>> G;
    vector<vector<ll>> flow;
    vector<int> prev_v, prev_e;
    
    MinCostFlow(int n): n(n) {
        G.assign(n, {});
        flow.assign(n, vector<ll>(n, 0));
        pot.assign(n, 0); dis.assign(n, 0);
        prev_v.assign(n, 0);prev_e.assign(n, 0);
    }
    void add_edge(int u, int v, ll cap, ll cost) {
        G[u].push_back({v, (int) G[v].size(), cap, cost});
        G[v].push_back({u, (int) G[u].size() - 1, 0, -cost});
    }
    pair<ll, ll> min_cost_flow(int S, int T, ll K) {
        ll flow_sent = 0, cost = 0;
        // If you have any negative-cost forward edges, uncomment this Bellman-Ford:
        // fill(pot.begin(), pot.end(), LLONG_MAX); pot[S] = 0;
        // for (int it = 0; it < n - 1; it++) for (int u = 0; u < n; u++)
        //     if (pot[u] < LLONG_MAX) for (auto &e : G[u])
        //         if (e.cap > 0) pot[e.to] = min(pot[e.to], pot[u] + e.cost);
        // otherwise, just zero them:        
        fill(pot.begin(), pot.end(), 0LL);
        while (flow_sent < K) {
            priority_queue<pair<ll, int>, vector<pair<ll, int>>, greater<pair<ll, int>>> pq;
            fill(dis.begin(), dis.end(), LLONG_MAX);
            dis[S] = 0; pq.push({0, S});
            while (pq.size()) {
                auto [d, u] = pq.top(); pq.pop();
                if (d > dis[u]) continue;
                for (int i = 0; i < (int) G[u].size(); i++) {
                    Edge &e = G[u][i];
                    if (e.cap > 0) {
                        ll nd = d + e.cost + pot[u] - pot[e.to];
                        if (nd < dis[e.to]) {
                            dis[e.to] = nd; prev_v[e.to] = u; prev_e[e.to] = i;
                            pq.push({nd, e.to});
                        }
                    }
                }
            }
            if (dis[T] == LLONG_MAX) break;
            for (int i = 0; i < n; i++) 
                if (dis[i] < LLONG_MAX) pot[i] += dis[i];
            ll addf = K - flow_sent;
            for (int v = T; v != S; v = prev_v[v]) 
                addf = min(addf, G[prev_v[v]][prev_e[v]].cap);
            for (int v = T; v != S; v = prev_v[v]) {
                auto &e = G[prev_v[v]][prev_e[v]];
                e.cap -= addf; G[v][e.rev].cap += addf;
                flow[prev_v[v]][v] += addf, flow[v][prev_v[v]] -= addf;
            }
            flow_sent += addf;
            cost += addf * pot[T];
        }
        return {flow_sent, cost};
    }
};

struct LowerBoundMCF {
    int n, SS, TT;
    MinCostFlow *mcf;
    vector<ll> demand;

    LowerBoundMCF(int n): n(n), SS(n), TT(n + 1) {
        demand.assign(n, 0);
        mcf = new MinCostFlow(n + 2);
    }
    void add_edge(int u, int v, ll low, ll up, ll cost) {
        demand[u] -= low;
        demand[v] += low;
        mcf->add_edge(u, v, up - low, cost);
    }
    tuple<bool, ll, ll> feasible() {
        ll total = 0;
        for (int i = 0; i < n; ++i) {
            if (demand[i] > 0) {
                mcf->add_edge(SS, i, demand[i], 0);
                total += demand[i];
            } else if (demand[i] < 0) mcf->add_edge(i, TT, -demand[i], 0);
        }
        auto [flow, cost] = mcf->min_cost_flow(SS, TT, total);
        return {flow == total, flow, cost};
    }
    pair<ll,ll> solve(int S, int T, ll K) {
        mcf->add_edge(T, S, LLONG_MAX / 4, 0);
        bool is_feasible; ll flow, cost;
        tie(is_feasible, flow, cost) = feasible();
        if (!is_feasible) return {-1, -1};   
        return {flow, cost};
    }
};

/*============================================================================================================
Max‑Flow Min‑Cut Theorem
  • An s–t cut is any partition of vertices into S (containing s) and T (containing t)
  • The capacity of the cut is the sum of capacities of all edges from S to T
  • Any s→t flow cannot exceed the capacity of any s–t cut ⇒ max flow ≤ cut capacity
  • The Max‑Flow Min‑Cut Theorem states:
    “The value of a maximum s→t flow equals the capacity of a minimum s–t cut”
  • After computing max flow with Ford‑Fulkerson, you can find a min cut by:
    1) Running DFS/BFS from s in the final residual graph following edges with residual > 0
    2) Let S be all reachable vertices, T the rest, (S,T) is a minimum cut
============================================================================================================*/
