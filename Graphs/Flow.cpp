#include <bits/stdc++.h>
using namespace std;
using ll = long long int;

/*============================================================================================================
Flow Algorithms Toolkit (Max-Flow, Min-Cut, Min-Cost, Lower Bounds)

Overview
  • This file collects the most important flow algorithms used in competitive programming
  • Covers:
    – Maximum flow (multiple classic algorithms)
    – Feasible flow with lower / upper bounds
    – Min-cost max-flow (with and without lower bounds)
  • All implementations correctly handle:
    – Parallel edges
    – Bidirectional edges
    – Self-loops
  • Written for correctness first, performance second (contest-safe)

1) Edmonds–Karp (Ford–Fulkerson + BFS):

Description:
  • Repeatedly finds the shortest (fewest-edge) augmenting path using BFS
  • Augments flow along that path until no augmenting path exists

Order:
  • O(V·E²) worst-case

When to use:
  • Very small graphs
  • Feasibility checking (especially with lower bounds)
  • Educational / debugging purposes

Important notes:
  • Parallel edges accumulate capacity
  • Reverse edges are maintained explicitly in residual graph
  • Self-loops are inert (never used in BFS)

2) Push–Relabel (FIFO Preflow–Push):

Description:
  • Maintains a preflow and fixes excess locally using push + relabel
  • Uses height labels instead of global augmenting paths

Order:
  • O(V²·E) worst-case
  • Very fast in practice on dense graphs

When to use:
  • Dense graphs
  • Large capacities
  • When Edmonds–Karp is too slow and Dinic struggles

Important notes:
  • Heights increase monotonically
  • Self-loops never satisfy height condition → ignored safely

3) Dinic’s Algorithm:

Description:
  • Builds a level graph using BFS
  • Finds blocking flows using DFS on level graph

Order:
  • Typical: O(E·√V)
  • Unit capacities / bipartite matching: O(√V·E)
  • Worst-case: O(V².E) (rare in practice)

When to use:
  • Default choice for max-flow in contests
  • Bipartite matching
  • Medium to large sparse graphs

Important notes:
  • Parallel edges handled naturally as separate edge objects
  • Self-loops cannot appear in level graph → harmless

4) Feasible Flow with Lower Bounds (Edmonds–Karp based):

Problem:
  • Each edge has:
    – lower bound L
    – upper bound C
  • Need any flow satisfying all constraints

Technique:
  • Convert lower bounds into node demands
  • Add super-source (SS) and super-sink (TT)
  • Add edge T → S with infinite capacity
  • Run max-flow from SS to TT

Order:
  • Same as Edmonds–Karp on (V + 2) nodes

When to use:
  • Circulation with demands
  • Scheduling with minimum quotas
  • Preprocessing for min-cost flow with bounds

Important notes:
  • If max-flow < total demand → no feasible solution exists
  • Flow recovery preserves original lower bounds

5) Min-Cost Max-Flow (Successive Shortest Path + Potentials):

Description:
  • Sends flow while minimizing total cost
  • Uses Dijkstra with Johnson potentials to avoid negative edges

Order:
  • O(F·E·log(V)), where F = amount of flow sent

When to use:
  • Cost-sensitive routing
  • Assignment / transportation problems
  • Minimum-cost circulation

Important notes:
  • Potentials ensure non-negative reduced costs
  • Bellman-Ford is only needed if negative cost edges exist initially
  • Self-loops do not affect shortest paths

6) Min-Cost Flow with Lower Bounds:

Description:
  • Combines:
    – Lower-bound feasibility
    – Cost minimization
  • Two-phase solution:
    1) Check feasibility via SS / TT
    2) Push additional flow from S → T at minimum cost

Order:
  • Feasibility: O(E·log(V))
  • Optimization: O(F·E·log(V))

When to use:
  • Scheduling with minimum + maximum limits
  • Resource allocation with costs
  • Advanced flow problems (hard CP tasks)

7) Max-Flow Min-Cut Theorem:

Core Definitions:
  • An s–t cut is any partition of the vertex set V into two disjoint sets (S, T) such that:
    – s ∈ S
    – t ∈ T
  • The capacity of a cut (S, T) is:
      sum of capacities of all directed edges u→v where u ∈ S and v ∈ T

Fundamental Properties:
  • Any s→t flow must cross from S to T somewhere
  • Therefore, the value of any feasible flow is upper-bounded by the capacity of every s–t cut:
      max_flow ≤ capacity(S, T)   for all cuts (S, T)

The Theorem (Strong Form):
  • The maximum possible value of an s→t flow is exactly equal to the minimum cut capacity:
      max_flow = min_cut

Residual-Graph Interpretation:
  • After running a max-flow algorithm:
    – Consider the final residual graph
    – Perform DFS/BFS from source s using only edges with residual capacity > 0
  • Let:
    – S = all vertices reachable from s
    – T = all remaining vertices
  • Then:
    – (S, T) forms a minimum s–t cut
    – All edges from S to T are saturated (residual = 0)
    – All edges from T to S have positive residual (flow can be canceled)

Takeaway:
  • Max-flow and min-cut are two sides of the same coin
  • Computing one automatically gives you the other
  • This duality is one of the most powerful ideas in graph algorithms
============================================================================================================*/

struct EdmondsKarp {
    int n;
    vector<vector<int>> G;
    vector<vector<ll>> cap, flow;

    explicit EdmondsKarp(const int n): n(n) {
        G.assign(n, {});
        cap.assign(n, vector<ll>(n, 0));
        flow.assign(n, vector<ll>(n, 0));
    }
    void add_edge(const int u, const int v, const ll c) {
        if (cap[u][v] == 0 and cap[v][u] == 0)
            G[u].push_back(v), G[v].push_back(u);
        cap[u][v] += c;
    }
    ll max_flow(int S, const int T, const ll K = LLONG_MAX) {
        ll total = 0; vector<int> par(n);
        while (total < K) {
            fill(par.begin(), par.end(), -1);
            par[S] = -2; queue<pair<int, ll>> q;
            q.emplace(S, LLONG_MAX); ll new_flow = 0;
            while (!q.empty()) {
                auto [u, f] = q.front(); q.pop();
                for (int v : G[u])
                    if (par[v] < 0 and cap[u][v] > 0) {
                        par[v] = u; ll mf = min(f, cap[u][v]);
                        if (v == T) { new_flow = mf; break; }
                        q.emplace(v, mf);
                    }
                if (new_flow) break;
            }
            if (!new_flow) break;
            new_flow = min(new_flow, K - total);
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

struct PushRelabel {
    int n;
    vector<ll> excess;
    queue<int> active;
    vector<int> height, seen;
    vector<vector<ll>> cap, flow;

    explicit PushRelabel(const int n): n(n) {
        excess.assign(n, 0);
        cap.assign(n, vector<ll>(n, 0));
        flow.assign(n, vector<ll>(n, 0));
        seen.assign(n, 0), height.assign(n, 0);
    }
    void add_edge(const int u, const int v, const ll c) {
        cap[u][v] += c;
    }
    void push(const int u, const int v, const int T) {
        const ll send = min(excess[u], cap[u][v] - flow[u][v]);
        if (send <= 0 or height[u] != height[v] + 1) return;
        excess[u] -= send, excess[v] += send;
        flow[u][v] += send, flow[v][u] -= send;
        if (excess[v] == send and v != T) active.push(v);
    }
    void relabel(const int u) {
        ll mh = LLONG_MAX;
        for (int v = 0; v < n; v++)
            if (cap[u][v] - flow[u][v] > 0)
                mh = min(mh, static_cast<ll>(height[v]));
        if (mh < LLONG_MAX) height[u] = static_cast<int>(mh) + 1;
    }
    void discharge(const int u, const int T) {
        while (excess[u] > 0)
            if (seen[u] < n) { const int v = seen[u]++; push(u, v, T); }
            else { relabel(u); seen[u] = 0; }
    }
    ll max_flow(const int S, const int T, const ll K = LLONG_MAX) {
        height.assign(n, 0);
        height[S] = n;
        ll flow_sent = 0;
        for (int u = 0; u < n; u++) {
            if (cap[S][u] > 0) {
                const ll want = min(cap[S][u], K - flow_sent);
                if (want <= 0) break;
                flow[S][u] = want; flow[u][S] = -want;
                excess[u] = want; excess[S] -= want;
                flow_sent += want;
                if (u != T) active.push(u);
                if (flow_sent == K) return flow_sent;
            }
        }
        while (!active.empty()) {
            const int u = active.front(); active.pop();
            if (u != S and u != T) discharge(u, T);
            if (excess[T] >= K) return K;
        }
        return min(excess[T], K);
    }
};

struct Dinic {
    struct Edge {
        int to, rev; ll cap;
    };

    int n, S, T;
    vector<int> level, ptr;
    vector<vector<Edge>> G;
    vector<vector<ll>> flow;
    vector<tuple<int, int, ll>> edges;

    Dinic(const int n, const int S, const int T): n(n), S(S), T(T) {
        flow.assign(n, vector<ll>(n, 0));
        G.assign(n, {}); level.assign(n, 0); ptr.assign(n, 0);
    }
    void add_edge(int u, int v, ll cap) {
        G[u].push_back({v, static_cast<int>(G[v].size()), cap});
        G[v].push_back({u, static_cast<int>(G[u].size()) - 1, 0});
        edges.emplace_back(u, v, cap);
    }
    bool BFS() {
        fill(level.begin(), level.end(), -1);
        level[S] = 0; queue<int> q; q.push(S);
        while (!q.empty()) {
            const int u = q.front(); q.pop();
            for (auto &e : G[u])
                if (level[e.to] < 0 and e.cap > 0) {
                    level[e.to] = level[u] + 1;
                    q.push(e.to);
                }
        }
        return level[T] >= 0;
    }
    ll DFS(const int u, const ll pushed) {
        if (!pushed or u == T) return pushed;
        for (int &cid = ptr[u]; cid < static_cast<int>(G[u].size()); cid++) {
            auto &[to, rev, cap] = G[u][cid];
            if (level[to] != level[u] + 1 or cap == 0) continue;
            if (const ll tr = DFS(to, min(pushed, cap))) {
                cap -= tr; G[to][rev].cap += tr;
                flow[u][to] += tr; flow[to][u] -= tr;
                return tr;
            }
        }
        return 0;
    }
    ll max_flow(const ll K = LLONG_MAX) {
        ll flow = 0;
        while (flow < K and BFS()) {
            fill(ptr.begin(), ptr.end(), 0);
            while (flow < K) {
                const ll pushed = DFS(S, K - flow);
                if (!pushed) break;
                flow += pushed;
            }
        }
        return flow;
    }
};

struct FeasibleFlow {
    int n;
    vector<vector<int>> G;
    vector<vector<ll>> cap, flow;

    explicit FeasibleFlow(const int n): n(n) {
        G.assign(n, {});
        cap.assign(n, vector<ll>(n, 0));
        flow.assign(n, vector<ll>(n, 0));
    }
    void add_edge(const int u, const int v, const ll c) {
        if (cap[u][v] == 0 and cap[v][u] == 0) {
            G[u].push_back(v);
            G[v].push_back(u);
        }
        cap[u][v] += c;
    }
    ll BFS(int S, const int T, vector<int> &par) const {
        fill(par.begin(), par.end(), -1);
        par[S] = -2; queue<pair<int, ll>> q;
        q.emplace(S, LLONG_MAX);
        while (!q.empty()) {
            auto [u, f] = q.front(); q.pop();
            for (int v : G[u])
                if (par[v] == -1 and cap[u][v] > 0) {
                    par[v] = u; ll newf = min(f, cap[u][v]);
                    if (v == T) return newf;
                    q.emplace(v, newf);
                }
        }
        return 0;
    }
    ll max_flow(const int S, const int T, const ll K = LLONG_MAX) {
        vector<int> par(n);
        ll flow_sum = 0, newf;
        while (flow_sum < K && (newf = BFS(S, T, par)) > 0) {
            newf = min(newf, K - flow_sum);
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
    [[nodiscard]] ll solve(const int m, const int S, const int T) const {
        vector<ll> balance(n, 0);
        vector<tuple<int,int,ll,ll>> edges;
        for (int i = 0; i < m; i++) {
            int u, v; ll L, C;
            cin >> u >> v >> L >> C;
            edges.emplace_back(--u, --v, L, C);
            balance[u] -= L; balance[v] += L;
        }
        const int SS = n, TT = n + 1; ll totalSupply = 0;
        FeasibleFlow ff(n + 2);
        for (int i = 0; i < n; i++)
            if (balance[i] > 0) {
                ff.add_edge(SS, i, balance[i]);
                totalSupply += balance[i];
            } else if (balance[i] < 0) ff.add_edge(i, TT, -balance[i]);
        for (const auto &e: edges) {
            int u, v; ll L, C; tie(u, v, L, C) = e;
            ff.add_edge(u, v, C - L);
        }
        ff.add_edge(T, S, LLONG_MAX);
        if (const ll fl = ff.max_flow(SS, TT); fl != totalSupply) {
            cout << "No feasible flow exists\n";
            return - 1;
        }
        return ff.max_flow(S, T);
    }
};

struct MinCostFlow {
    struct Edge {
        int to, rev; ll cap, cost;
    };

    int n;
    vector<ll> dis, pot;
    vector<vector<Edge>> G;
    vector<vector<ll>> flow;
    vector<int> prev_v, prev_e;

    explicit MinCostFlow(const int n): n(n) {
        G.assign(n, {});
        flow.assign(n, vector<ll>(n, 0));
        pot.assign(n, 0); dis.assign(n, 0);
        prev_v.assign(n, 0);prev_e.assign(n, 0);
    }
    void add_edge(const int u, const int v, const ll cap, const ll cost) {
        G[u].push_back({v, static_cast<int>(G[v].size()), cap, cost});
        G[v].push_back({u, static_cast<int>(G[u].size()) - 1, 0, -cost});
    }
    pair<ll, ll> min_cost_flow(int S, const int T, const ll K = LLONG_MAX) {
        ll flow_sent = 0, cost = 0;
        // If you have any negative-cost forward edges, uncomment this Bellman-Ford:
        // fill(pot.begin(), pot.end(), LLONG_MAX); pot[S] = 0;
        // for (int it = 0; it < n - 1; it++) for (int u = 0; u < n; u++)
        //     if (pot[u] < LLONG_MAX) for (auto &e : G[u])
        //         if (e.cap > 0) pot[e.to] = min(pot[e.to], pot[u] + e.cost);
        // otherwise, just zero them:
        fill(pot.begin(), pot.end(), 0LL);
        while (flow_sent < K) {
            priority_queue<pair<ll, int>, vector<pair<ll, int>>, greater<>> pq;
            fill(dis.begin(), dis.end(), LLONG_MAX);
            dis[S] = 0; pq.emplace(0, S);
            while (!pq.empty()) {
                auto [d, u] = pq.top(); pq.pop();
                if (d > dis[u]) continue;
                for (int i = 0; i < static_cast<int>(G[u].size()); i++) {
                    if (Edge &e = G[u][i]; e.cap > 0) {
                        if (ll nd = d + e.cost + pot[u] - pot[e.to]; nd < dis[e.to]) {
                            dis[e.to] = nd; prev_v[e.to] = u; prev_e[e.to] = i;
                            pq.emplace(nd, e.to);
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

    explicit LowerBoundMCF(const int n): n(n), SS(n), TT(n + 1) {
        demand.assign(n, 0);
        mcf = new MinCostFlow(n + 2);
    }
    void add_edge(const int u, const int v, const ll low, const ll up, const ll cost) {
        demand[u] -= low;
        demand[v] += low;
        mcf->add_edge(u, v, up - low, cost);
    }
    [[nodiscard]] tuple<bool, ll, ll> feasible() const {
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
    [[nodiscard]] pair<ll,ll> solve(const int S, const int T, ll K = LLONG_MAX) const {
        mcf->add_edge(T, S, LLONG_MAX / 4, 0);
        auto [is_feasible, flow, cost] = feasible();
        if (!is_feasible) return {-1, -1};
        return {flow, cost};
    }
};
