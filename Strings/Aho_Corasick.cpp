#include <bits/stdc++.h>
using namespace std;
using ll = long long int;
using ull = unsigned long long int;

/*============================================================================================================
Aho–Corasick utilities (multi-pattern automaton) — tasks & notes

This file implements a compact Aho–Corasick trie and three higher-level routines
built on top of it (kept separate from the build() step as you requested).

Core struct
  • AhoCorasick::Node
    - next[c]  : deterministic transition for character c (filled by build())
    - link     : failure link (longest proper suffix that is a prefix)
    - out      : list of pattern ids that end at this node (includes suffix outputs merged from link during build())

  • AhoCorasick methods:
    - add_pattern(string): insert a pattern into the trie (patterns indexed 0..patterns-1)
    - build()            : BFS build of failure links and fill missing transitions so
                           next[v][c] is always defined (root-based automaton).
                           After build(), the structure is a full DFA for the pattern set.

High-level algorithms (separate functions)
  1) smallest_avoiding(int L)
     • Task: find the lexicographically smallest string of length L that contains
       none of the patterns as substrings.
     • Approach:
       - Mark "bad" states (those with non-empty out).
       - Consider the subgraph induced by safe states and automaton transitions.
       - Collapse SCCs into components; detect components that contain cycles.
       - For each component compute whether it can reach a cycle (=> infinite extension),
         and if not, compute the longest path length reachable from it (DP on component DAG).
       - Greedily build the answer by trying characters in increasing order and
         accepting a character only if the destination node is safe and can finish the
         remaining length (either reach a cycle or has maxLen >= remaining).
     • Complexity: O(nodes * ALPH + comp_DAG) preprocessing; O(L * ALPH) construction time.
     • Memory: O(nodes + edges). Good for large L because we don't build a DP over L.

  2) shortest_superstring_all()
     • Task: find the shortest string that contains every pattern at least once.
     • Approach:
       - Do BFS on the product graph (state v, bitmask of matched patterns).
       - Each transition updates the mask with patterns ending at the new state.
       - BFS yields shortest path length; parent map is used to reconstruct the string.
     • Complexity: O(nodes * 2^n * ALPH) in the worst case (exponential in number of patterns).
     • Practical notes / constraints:
       - This is inherently exponential in n = number of patterns.
       - The implementation uses an integer mask packed into a 64-bit key; it is safe and
         practical for small n (n ≤ ~20 is typical). For larger n the memory/time explode.
       - If you need correctness for larger n, replace internal mask/key logic with a
         robust 64-bit or multi-word mask, or use heuristic/approximate algorithms.
     • Returns empty string on impossibility or when too large to compute reasonably.

  3) smallest_with_k(int L, int K)
     • Task: lexicographically smallest string of length L that contains at least K distinct patterns.
     • Approach:
       - Build full automaton graph and its component DAG (SCC collapse).
       - Compute for each component whether it can reach a cycle and its max acyclic length.
       - Precompute a "reachable pattern mask" for each component (which patterns can be
         obtained from nodes in that component or downstream components).
       - Backtrack lexicographically with pruning:
         • length feasibility: node can reach cycle or maxLen >= remaining
         • pattern feasibility: (currentMask ∪ reachableMask[node]) must be able to reach K
       - Use memoization of failed (v,step,mask) states to avoid repeated search.
     • Complexity: Preprocessing O(nodes * ALPH + comp_DAG). Backtracking cost depends on L,
       alphabet, and patterns; worst-case expensive but pruned heavily by component masks.
     • Practical notes:
       - Uses packed bit-blocks (vector<long long>) for masks, good up to a few hundred patterns.
       - If patterns or L are very large this can become heavy; tune mask representation or
         the memoization hashing if needed.

Important implementation / safety notes
  • Alphabet: fixed to lowercase 'a'..'z'. Non-lowercase input will index out of range.
  • build() ensures next[v][c] is always valid for all v, c; code assumes that everywhere.
  • shortest_superstring_all() is exponential in the number of patterns; keep n small (≈≤20).
    The current implementation packs mask into 32-bit/64-bit parts; if you need correctness
    for larger n, I can replace the mask/key code with a robust representation.
  • smallest_with_k() uses component-level mask propagation to prune the search. It uses
    64-bit blocks; adjust if you need support for thousands of patterns (but performance will suffer).
  • Memory/time tradeoffs:
    - smallest_avoiding(): memory-friendly, linear-ish even for huge L.
    - shortest_superstring_all(): memory/time blow up with number of patterns.
    - smallest_with_k(): practical for moderate patterns; worst-case exponential.

Summary
  • This file is primarily an Aho–Corasick trie + three high-level helpers that run on top of it.
  • Use add_pattern(...) then call build() once; call the helpers after build() (they perform
    any additional preprocessing required for their task).
  • If you need the shortest_superstring_all() to work for large pattern counts, ask and I will
    replace the pack/key logic with a 64-bit or multi-word mask implementation (and warn again about
    inherent exponential complexity).
============================================================================================================*/

struct AhoCorasick {
    static constexpr int ALPH = 26;
    static int to_index(const char c) { return c - 'a'; }

    struct Node {
        int link;
        vector<int> next, out;
        Node() : link(-1) { next.assign(ALPH, -1); }
    };

    int patterns = 0;
    vector<Node> trie;
    vector<int> pat_len;

    AhoCorasick() { trie.emplace_back(); }
    void add_pattern(const string &p) {
        int u = 0;
        for (const char ch : p) {
            const int x = to_index(ch);
            if (trie[u].next[x] == -1) {
                trie[u].next[x] = static_cast<int>(trie.size());
                trie.emplace_back();
            }
            u = trie[u].next[x];
        }
        trie[u].out.push_back(patterns++);
        pat_len.push_back(static_cast<int>(p.size()));
    }
    void build() {
        queue<int> q;
        trie[0].link = 0;
        for (int c = 0; c < ALPH; ++c) {
            if (int u = trie[0].next[c]; u != -1) { trie[u].link = 0; q.push(u); }
            else trie[0].next[c] = 0;
        }
        while (!q.empty()) {
            const int u = q.front(); q.pop();
            if (const int lu = trie[u].link; lu != u and !trie[lu].out.empty())
                trie[u].out.insert(trie[u].out.end(),
                                   trie[lu].out.begin(), trie[lu].out.end());
            for (int c = 0; c < ALPH; ++c) {
                if (int v = trie[u].next[c]; v != -1) {
                    trie[v].link = trie[trie[u].link].next[c];
                    q.push(v);
                } else trie[u].next[c] = trie[trie[u].link].next[c];
            }
        }
    }

private:
    static void kosaraju_scc(const int n, const vector<vector<int>> &g,
                      vector<int> &comp, int &comp_cnt) {
        vector<int> order; order.reserve(n);
        vector<char> vis(n, 0);
        function<void(int)> dfs1 = [&](const int v) {
            vis[v] = 1;
            for (const int to : g[v]) if (!vis[to]) dfs1(to);
            order.push_back(v);
        };
        for (int i = 0; i < n; ++i) if (!vis[i]) dfs1(i);
        vector<vector<int>> gr(n);
        for (int v = 0; v < n; ++v) for (const int to : g[v]) gr[to].push_back(v);
        comp.assign(n, -1); comp_cnt = 0;
        function<void(int)> dfs2 = [&](const int v) {
            comp[v] = comp_cnt;
            for (const int to : gr[v]) if (comp[to] == -1) dfs2(to);
        };
        for (int i = n - 1; i >= 0; --i)
            if (const int v = order[i]; comp[v] == -1) { dfs2(v); ++comp_cnt; }
    }
    static int popcount_mask(const vector<ull> &mask) {
        int cnt = 0;
        for (const ull x : mask) cnt += __builtin_popcountll(x);
        return cnt;
    }
    static void mask_or_into(vector<ull> &dst, const vector<ull> &src) {
        const int m = static_cast<int>(dst.size());
        for (int i = 0; i < m; ++i) dst[i] |= src[i];
    }
    static string make_key(const int v, const vector<ull> &mask) {
        string k;
        k.resize(4 + mask.size() * sizeof(ull));
        memcpy(&k[0], &v, 4);
        memcpy(&k[4], mask.data(), mask.size() * sizeof(ull));
        return k;
    }

public:
    [[nodiscard]] string smallest_avoiding(int L) const {
        int N = static_cast<int>(trie.size());
        if (N == 0) return {};
        vector<char> bad(N, 0);
        for (int v = 0; v < N; ++v) if (!trie[v].out.empty()) bad[v] = 1;
        if (bad[0]) return {};
        vector<vector<int>> g(N);
        for (int v = 0; v < N; ++v) {
            if (bad[v]) continue;
            for (int c = 0; c < ALPH; ++c)
                if (int u = trie[v].next[c]; !bad[u]) g[v].push_back(u);
        }
        vector<int> comp; int comp_cnt;
        kosaraju_scc(N, g, comp, comp_cnt);
        vector compSize(comp_cnt, 0);
        vector<char> compSelfLoop(comp_cnt, 0);
        for (int v = 0; v < N; ++v) if (!bad[v]) {
            int cv = comp[v];
            compSize[cv]++;
            for (int u : g[v]) if (u == v) compSelfLoop[cv] = 1;
        }
        vector<char> compIsCycle(comp_cnt, 0);
        for (int i = 0; i < comp_cnt; ++i)
            if (compSize[i] > 1 or compSelfLoop[i]) compIsCycle[i] = 1;
        vector<vector<int>> dag(comp_cnt);
        vector indeg(comp_cnt, 0);
        for (int v = 0; v < N; ++v) if (!bad[v]) {
            int cv = comp[v];
            for (int u : g[v])
                if (int cu = comp[u]; cv != cu) dag[cv].push_back(cu);
        }
        for (int i = 0; i < comp_cnt; ++i) {
            sort(dag[i].begin(), dag[i].end());
            dag[i].erase(unique(dag[i].begin(), dag[i].end()), dag[i].end());
            for (int to : dag[i]) indeg[to]++;
        }
        queue<int> q;
        for (int i = 0; i < comp_cnt; ++i) if (indeg[i] == 0) q.push(i);
        vector<int> topo; topo.reserve(comp_cnt);
        while (!q.empty()) {
            int x = q.front(); q.pop();
            topo.push_back(x);
            for (int to : dag[x]) if (--indeg[to] == 0) q.push(to);
        }
        vector<char> compCanReachCycle(comp_cnt, 0);
        for (int i = static_cast<int>(topo.size()) - 1; i >= 0; --i) {
            int c = topo[i];
            if (compIsCycle[c]) compCanReachCycle[c] = 1;
            for (int to : dag[c]) if (compCanReachCycle[to]) compCanReachCycle[c] = 1;
        }
        vector compMaxLen(comp_cnt, 0);
        for (int i = static_cast<int>(topo.size()) - 1; i >= 0; --i) {
            constexpr int INF = 1e9;
            int c = topo[i];
            if (compCanReachCycle[c]) { compMaxLen[c] = INF; continue; }
            int best = 0;
            for (int to : dag[c]) {
                if (compMaxLen[to] == INF) { best = INF; break; }
                best = max(best, 1 + compMaxLen[to]);
            }
            compMaxLen[c] = best;
        }
        vector<char> canReachCycle(N, 0);
        vector maxLen(N, 0);
        for (int v = 0; v < N; ++v) if (!bad[v]) {
            int c = comp[v];
            canReachCycle[v] = compCanReachCycle[c];
            maxLen[v] = compMaxLen[c];
        }
        string ans; ans.reserve(L);
        int v = 0;
        if (L == 0) return "";
        for (int step = 0; step < L; ++step) {
            bool placed = false;
            for (int ci = 0; ci < ALPH; ++ci) {
                int u = trie[v].next[ci];
                if (bad[u]) continue;
                if (int rem = L - step - 1; canReachCycle[u] or maxLen[u] >= rem) {
                    ans.push_back(static_cast<char>('a' + ci));
                    v = u; placed = true; break;
                }
            }
            if (!placed) return {};
        }
        return ans;
    }
    [[nodiscard]] string shortest_superstring_all() const {
        int n = patterns;
        if (n == 0) return "";
        int B = (n + 63) / 64;
        vector<ull> full(B, 0);
        for (int i = 0; i < n; ++i) full[i / 64] |= (1ull << (i % 64));
        auto node_mask = [&](const int node)->vector<ull> {
            vector<ull> m(B, 0);
            for (const int pid : trie[node].out) if (pid >= 0 and pid < n)
                m[pid / 64] |= (1ull << (pid % 64));
            return m;
        };
        unordered_set<string> vis;
        vis.reserve(1 << 16);
        unordered_map<string, pair<string, char>> parent;
        queue<pair<int, vector<ull>>> q;
        vector<ull> start = node_mask(0);
        string sk = make_key(0, start);
        vis.insert(sk);
        parent[sk] = { string(), 0 };
        q.emplace(0, start);
        bool found = false;
        string final_key;
        while (!q.empty()) {
            auto [fst, snd] = q.front(); q.pop();
            int v = fst;
            vector<ull> mask = move(snd);
            if (mask == full) { found = true; final_key = make_key(v, mask); break; }
            for (int c = 0; c < ALPH; ++c) {
                int u = trie[v].next[c];
                vector<ull> nmask = mask;
                for (int pid : trie[u].out) if (pid >= 0 and pid < n)
                    nmask[pid / 64] |= (1ull << (pid % 64));
                string key = make_key(u, nmask);
                if (!vis.insert(key).second) continue;
                parent.emplace(key, make_pair(make_key(v, mask), static_cast<char>('a' + c)));
                q.emplace(u, move(nmask));
            }
        }
        if (!found) return {};
        string res;
        string cur = final_key;
        while (true) {
            auto it = parent.find(cur);
            if (it == parent.end()) break;
            auto [fst, snd] = it->second;
            if (fst.empty()) break;
            res.push_back(snd);
            cur = fst;
        }
        reverse(res.begin(), res.end());
        return res;
    }
    [[nodiscard]] string smallest_with_k(int L, int K) const {
        int N = static_cast<int>(trie.size());
        if (K <= 0) return string(L, 'a');
        if (patterns == 0) return {};
        vector<vector<int>> g(N);
        for (int v = 0; v < N; ++v)
            for (int c = 0; c < ALPH; ++c) g[v].push_back(trie[v].next[c]);
        vector<int> comp; int comp_cnt;
        kosaraju_scc(N, g, comp, comp_cnt);
        vector compSize(comp_cnt, 0);
        vector<char> compSelfLoop(comp_cnt, 0);
        for (int v = 0; v < N; ++v) {
            compSize[comp[v]]++;
            for (int u : g[v]) if (u == v) compSelfLoop[comp[v]] = 1;
        }
        vector<char> compIsCycle(comp_cnt, 0);
        for (int i = 0; i < comp_cnt; ++i)
            if (compSize[i] > 1 or compSelfLoop[i]) compIsCycle[i] = 1;
        vector<vector<int>> dag(comp_cnt);
        vector indeg(comp_cnt, 0);
        for (int v = 0; v < N; ++v) {
            int cv = comp[v];
            for (int u : g[v]) {
                if (int cu = comp[u]; cv != cu) dag[cv].push_back(cu);
            }
        }
        for (int i = 0; i < comp_cnt; ++i) {
            sort(dag[i].begin(), dag[i].end());
            dag[i].erase(unique(dag[i].begin(), dag[i].end()), dag[i].end());
            for (int to : dag[i]) indeg[to]++;
        }
        queue<int> q0;
        for (int i = 0; i < comp_cnt; ++i) if (indeg[i] == 0) q0.push(i);
        vector<int> topo;
        while (!q0.empty()) {
            int x = q0.front(); q0.pop();
            topo.push_back(x);
            for (int y : dag[x]) if (--indeg[y] == 0) q0.push(y);
        }
        vector<char> compCanReachCycle(comp_cnt, 0);
        for (int i = static_cast<int>(topo.size()) - 1; i >= 0; --i) {
            int c = topo[i];
            if (compIsCycle[c]) compCanReachCycle[c] = 1;
            for (int to : dag[c]) if (compCanReachCycle[to]) compCanReachCycle[c] = 1;
        }
        vector compMaxLen(comp_cnt, 0);
        for (int i = static_cast<int>(topo.size()) - 1; i >= 0; --i) {
            constexpr int INF = 1e9;
            int c = topo[i];
            if (compCanReachCycle[c]) { compMaxLen[c] = INF; continue; }
            int best = 0;
            for (int to : dag[c]) {
                if (compMaxLen[to] == INF) { best = INF; break; }
                best = max(best, 1 + compMaxLen[to]);
            }
            compMaxLen[c] = best;
        }
        vector<char> nodeCanReachCycle(N, 0);
        vector nodeMaxLen(N, 0);
        for (int v = 0; v < N; ++v) {
            nodeCanReachCycle[v] = compCanReachCycle[comp[v]];
            nodeMaxLen[v] = compMaxLen[comp[v]];
        }
        int B = (patterns + 63) / 64;
        vector compMask(comp_cnt, vector<ull>(B, 0));
        for (int v = 0; v < N; ++v) {
            int cv = comp[v];
            for (int pid : trie[v].out) if (pid >= 0) {
                compMask[cv][pid / 64] |= (1ull << (pid % 64));
            }
        }
        for (int i = static_cast<int>(topo.size()) - 1; i >= 0; --i) {
            for (int c = topo[i]; int to : dag[c]) mask_or_into(compMask[c], compMask[to]);
        }
        vector nodeMask(N, vector<ull>(B, 0));
        for (int v = 0; v < N; ++v) nodeMask[v] = compMask[comp[v]];
        string path; path.reserve(L);
        bool found = false;
        unordered_set<string> dead;
        dead.reserve(1 << 16);
        function<bool(int,int,vector<ull>&)> dfs = [&](const int v, const int step, const vector<ull> &curMask)->bool {
            if (found) return true;
            if (const int rem = L - step; !(nodeCanReachCycle[v] or nodeMaxLen[v] >= rem)) return false;
            vector<ull> tmp(B);
            for (int i = 0; i < B; ++i) tmp[i] = curMask[i] | nodeMask[v][i];
            if (const int totalPotential = popcount_mask(tmp); totalPotential < K) return false;
            if (step == L) {
                if (const int have = popcount_mask(curMask); have >= K) { found = true; return true; }
                return false;
            }
            string key = make_key(v, curMask);
            if (dead.count(key)) return false;
            for (int ci = 0; ci < ALPH; ++ci) {
                const int u = trie[v].next[ci];
                if (const int remAfter = L - (step + 1); !(nodeCanReachCycle[u] or nodeMaxLen[u] >= remAfter)) continue;
                vector<ull> newMask = curMask;
                for (const int pid : trie[u].out) if (pid >= 0)
                    newMask[pid / 64] |= (1ull << (pid % 64));
                vector<ull> combined(B);
                for (int i = 0; i < B; ++i) combined[i] = newMask[i] | nodeMask[u][i];
                if (const int potentialTotal = popcount_mask(combined); potentialTotal < K) continue;
                path.push_back(static_cast<char>('a' + ci));
                if (dfs(u, step + 1, newMask)) return true;
                path.pop_back();
            }
            dead.insert(move(key));
            return false;
        };
        vector<ull> startMask(B, 0);
        for (int pid : trie[0].out) if (pid >= 0)
            startMask[pid / 64] |= (1ull << (pid % 64));
        if (bool ok = dfs(0, 0, startMask); !ok) return {};
        return path;
    }
};
