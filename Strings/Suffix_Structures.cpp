#include <bits/stdc++.h>
#include <utility>
using namespace std;
using ll = long long int;

/*============================================================================================================
Suffix Data Structures

This file contains two core string data structures used in competitive programming and string algorithms:
  1) SuffixArray (classic doubling + Kasai LCP + RMQ)
  2) SuffixAutomaton (online SAM with clones)

1) SuffixArray:

Purpose
  • Provide lexicographic ordering of all suffixes of a string (suffix array).
  • Support pattern search (binary search on SA), LCP queries, and counting distinct substrings.
  • Provide two ways to get LCP between suffixes:
    a) Kasai LCP array + RMQ (constant-time LCP queries after preprocessing)
    b) LCP by comparing stored doubling-level class arrays (logarithmic time, no RMQ)

How it works (brief)
  • sort_cyclic_shifts: standard doubling algorithm (sort by 2^k-length prefixes using counting sort by classes)
    - Optionally stores equivalence classes at each doubling level (c_levels) so LCP can be computed by levels.
  • Constructor flow:
      1) Build class levels for the base string (for lcp_by_levels usage)
      2) Build full suffix array on s + sentinel (char(0)), then strip sentinel to get SA for original string
      3) Build inverse array pos (pos[i] = SA rank of suffix starting at i)
      4) Run Kasai to compute LCP array
      5) Build sparse-table RMQ over LCP for O(1) min queries
  • find_pattern uses binary search on SA and compare_suffix_pattern for lexicographic compare.
  • lcp_by_levels uses stored class arrays to quickly check equality of 2^k blocks.
  • lcp_query_by_lcp_array uses RMQ on Kasai LCP array to return exact LCP between two suffixes.

Supported operations (interface)
  • find_pattern(pat): returns vector of starting indices where pattern occurs (sorted)
  • lcp_by_levels(i,j): LCP using c_levels (O(log(n)))
  • lcp_query_by_lcp_array(i,j): LCP via pos + RMQ on Kasai LCP (O(1))
  • count_different_substrings(): total number of distinct substrings (using SA+LCP)

Complexities
  • sort_cyclic_shifts (doubling): O(n.log(n)) time, O(n) extra memory per level (radix counting)
  • Kasai LCP: O(n)
  • RMQ build (sparse table): O(n.log(n)) preprocessing, O(1) per query
  • find_pattern: O(|pat|.log(n)) comparisons across binary search steps
  • lcp_by_levels: O(log(n))
  • Memory:
    - SA, pos, LCP: O(n)
    - RMQ sparse table: O(n.log(n))
    - c_levels: O(n.log(n)) if stored for all levels

Important implementation notes
  • Uses sentinel char(0) to ensure proper ordering when building SA on s + sentinel.
  • sort_cyclic_shifts accepts a store flag to optionally save class arrays truncated to original n.
  • compare_suffix_pattern compares only up to pattern length and suffix end.
  • rmq returns 0 for empty ranges (l > r) to simplify callers.
  • All arrays are 0-indexed, suffix_array holds starting positions within [0, n-1].
  • Edge-cases: empty strings and single-character strings are naturally supported but watch calls using indices.

2) SuffixAutomaton:

Purpose (concise)
  - The suffix automaton (SAM) is the minimal deterministic automaton that recognizes all substrings of a given string S.
  - It compactly represents substring structure and supports membership, counting, enumeration, and many substring queries efficiently.

Core concepts
  • Language recognized: { x | x is a substring of S }.
  • Minimal DFA: SAM is the minimal deterministic finite automaton for that language (i.e., no smaller DFA exists).
  • Alternative names: DAWG (Directed Acyclic Word Graph), factor automaton (context dependent).

State semantics (the most important facts about states)
  • Endpos set:
    - Each state v corresponds to an equivalence class of substrings that share the same set of end-positions in S.
    - endpos(x) = { i | S[i - |x| + 1 .. i] == x }.
    - If two substrings have identical endpos sets they map to the same state.
  • Length interval:
    - All strings represented by state v have lengths in the interval (len[link[v]], len[v]].
    - len[v] = maximum length of strings in v.
    - minLength(v) = len[link[v]] + 1 (shortest length represented by v).
  • Transitions and right-contexts:
    - Outgoing labeled transitions from v correspond to characters that can follow substrings in v.
    - The set of outgoing labels is exactly the set of characters that appear as the immediate right-context after occurrences of any substring in v.
  • Suffix link:
    - link[v] points to the state of the longest proper suffix class of any string in v.
    - The suffix links form a rooted tree (link tree) rooted at state 0.
  • first_pos:
    - For many implementations, state.v.first_pos stores an index i meaning the longest string represented by v ends at i in S.
    - For a substring x at state v, a valid occurrence can be reconstructed as first_pos - |x| + 1.
  • cnt (occurrence count):
    - When constructed, newly created non-clone states get cnt = 1 (one new end position).
    - Cloned states get cnt = 0.
    - After construction, running a suffix-link DP (propagating counts from larger len to smaller len) aggregates real occurrence counts: cnt[v] = size(endpos-set).
  • Clones:
    - Clones are copies of an existing state q created when we need to split q because a new transition would violate len invariants.
    - clone.len = len[p] + 1 (where p is the state that caused the split).
    - clone receives a copy of q's transitions and link; clone.cnt = 0 (unless you later add counts).
    - q and cur have their suffix links redirected to clone so that length intervals and endpos semantics are preserved.
  • Uniqueness and mapping:
    - Each distinct substring x corresponds to exactly one pair (v, L) where v is a state and L is length such that len[link[v]] < L ≤ len[v].
    - Therefore the number of distinct substrings = Σ_{v != root} (len[v] - len[link[v]]).

Construction and invariants
  • Online extend(c) algorithm:
    1) Create new state cur with len = len[last] + 1 and cnt = 1.
    2) Walk back p = last via suffix links adding transition p->c = cur while missing.
    3) If p == -1: set link[cur] = 0.
    4) Else let q = p->c:
       - If len[p] + 1 == len[q]: set link[cur] = q.
       - Else create clone copying q's transitions and link; set clone.len = len[p] + 1; redirect links from q and cur to clone; fix transitions from p backwards that pointed to q to point to clone.
  • Guarantees:
    - Number of states ≤ 2n - 1.
    - Number of transitions is O(states * alphabet) but in practice O(n) overall (amortized).
    - len strictly increases along transitions (graph of transitions is acyclic in terms of len).
    - Link graph is a tree (each node has exactly one link parent except root).

Useful derived properties & formulas
  • Distinct substrings:
    - count = Σ_{v != 0} (len[v] - len[link[v]]).
  • Distinct substrings by length L:
    - For each state v, it contributes 1 to every length L in [len[link[v]]+1, len[v]].
    - So you can use a difference array across those intervals to compute counts per length in total O(states).
  • Occurrence counting:
    - After construction, process states in decreasing len order and add cnt[v] to cnt[link[v]].
    - Then cnt[v] equals number of end positions (occurrences) of any representative substring of v.
  • First/last positions:
    - Maintaining first_pos[u] when a state is created gives a quick way to compute some occurrence positions. last_pos can be similarly maintained by updating on insertion.
  • LCS (longest common substring) with another string T:
    - Walk T over the SAM: maintain current state u and match length l. If char missing, follow suffix links until possible. Track maximum l encountered.
    - Complexity O(|T|).
  • K-th lexicographic substring:
    - Build sorted outgoing transitions per state and dp[u] = Σ(1 + dp[v]) over outgoing edges u->v.
    - Walk from root in lexicographic order using dp blocks to pick k-th substring.
    - dp[0] equals total number of distinct substrings.
    - dp[v] = the number of distinct substrings whose path starts from state v.
  • Shortest absent word (minimal non-appearing string over an alphabet):
    - DP on states computing minimal missing length d[v]; reconstruct greedily.
  • Relationship with suffix tree:
    - SAM is closely related to suffix tree/DAWG; SAM compresses some suffix-tree structure into a minimal automaton. There's a mapping between tree edges and SAM states/paths (useful in proofs).
  • Reversed string:
    - Building SAM over reversed string is a convenient trick to handle prefix-based queries similarly (e.g., distinct prefixes, distinct substrings beginning at positions).

Common queries and required precomputation
  • contains(p): no extra precomputation; just walk transitions O(|p|).
  • occurrences_count(p): requires prepare_occurrence_info (aggregate cnt via suffix links) => O(|p|) to locate state, O(1) to return cnt.
  • first_occurrence(p): requires first_pos to be maintained (done at construction); return first_pos - |p| + 1.
  • all_occurrences(p): uses inv_link lists built after prepare_occurrence_info and DFS in suffix-link tree; deduplicate if needed.
  • kth_substring: requires prepare_kth (sorted transitions + dp).
  • LCS with T: no precomputation beyond SAM build, O(|T|).

When to use SAM vs SA
  • Use SAM when you need:
    - fast occurrence counting, first/last occurrences, enumeration (k-th), LCS with multiple strings, or an online build.
  • Use SA when you need:
    - suffix order, cheap LCP queries across many pairs, text indexing with small memory but heavier preprocessing.
============================================================================================================*/

struct SuffixArray {
    int n;
    string s;
    vector<vector<int>> c_levels, st;
    vector<int> suffix_array, pos, lcp, lg;

    explicit SuffixArray(string in): s(std::move(in)) {
        n = static_cast<int>(s.size());
        sort_cyclic_shifts(s, true);
        const string ss = s + static_cast<char>(0);
        suffix_array = sort_cyclic_shifts(ss, false);
        suffix_array.erase(suffix_array.begin());
        pos.assign(n, 0); for (int i = 0; i < n; i++) pos[suffix_array[i]] = i;
        build_kasai();
        build_rmq();
    }
    [[nodiscard]] vector<int> find_pattern(const string &pat) const {
        vector<int> res;
        auto cmp_at = [&](const int idx)->int {
            return compare_suffix_pattern(suffix_array[idx], pat);
        };
        int L = 0, R = n;
        while (L < R)
            if (const int mid = (L + R) >> 1; cmp_at(mid) < 0) L = mid + 1;
            else R = mid;
        const int lb = L;
        if (lb == n or cmp_at(lb) != 0) return res;
        L = 0; R = n;
        while (L < R)
            if (const int mid = (L + R) >> 1; cmp_at(mid) <= 0) L = mid + 1;
            else R = mid;
        const int ub = L;
        res.reserve(ub - lb);
        for (int i = lb; i < ub; i++) res.push_back(suffix_array[i]);
        sort(res.begin(), res.end());
        return res;
    }
    [[nodiscard]] int lcp_by_levels(int i, int j) const {
        if (i == j) return n - i;
        int ans = 0;
        for (int k = static_cast<int>(c_levels.size()) - 1; k >= 0; --k) {
            if (i >= n or j >= n) break;
            if (c_levels[k][i] == c_levels[k][j]) {
                const int step = 1 << k;
                ans += step;
                i += step; j += step;
                if (i > n or j > n) break;
            }
        }
        return ans;
    }
    [[nodiscard]] int lcp_query_by_lcp_array(const int i, const int j) const {
        if (i == j) return n - i;
        int pi = pos[i], pj = pos[j];
        if (pi > pj) swap(pi, pj);
        return rmq(pi, pj - 1);
    }
    [[nodiscard]] ll count_different_substrings() const {
        ll res = 0;
        for (int i = 0; i < n; i++) {
            res += (n - suffix_array[i]);
            if (i > 0) res -= lcp[i - 1];
        }
        return res;
    }

private:
    vector<int> sort_cyclic_shifts(const string &sref, const bool store) {
        const int m = static_cast<int>(sref.size()), alphabet = 256;
        vector<int> p(m), c(m), cnt(max(alphabet, m), 0);
        for (int i = 0; i < m; i++) cnt[static_cast<unsigned char>(sref[i])]++;
        for (int i = 1; i < alphabet; i++) cnt[i] += cnt[i - 1];
        for (int i = 0; i < m; i++) p[--cnt[static_cast<unsigned char>(sref[i])]] = i;
        c[p[0]] = 0; int classes = 1;
        for (int i = 1; i < m; ++i) {
            if (sref[p[i]] != sref[p[i - 1]]) ++classes;
            c[p[i]] = classes - 1;
        }
        vector<vector<int>> tmp_levels;
        if (store) tmp_levels.emplace_back(c.begin(), c.begin() + min(static_cast<int>(c.size()), n));
        vector<int> pn(m), cn(m);
        for (int h = 0; (1 << h) < m; h++) {
            for (int i = 0; i < m; i++) {
                pn[i] = p[i] - (1 << h);
                if (pn[i] < 0) pn[i] += m;
            }
            fill_n(cnt.begin(), classes, 0);
            for (int i = 0; i < m; i++) cnt[c[pn[i]]]++;
            for (int i = 1; i < classes; i++) cnt[i] += cnt[i - 1];
            for (int i = m - 1; i >= 0; i--) p[--cnt[c[pn[i]]]] = pn[i];
            cn[p[0]] = 0; classes = 1;
            for (int i = 1; i < m; i++) {
                pair cur = {c[p[i]], c[(p[i] + (1 << h)) % m]};
                if (pair prev = {c[p[i-1]], c[(p[i-1] + (1 << h)) % m]}; cur != prev) ++classes;
                cn[p[i]] = classes - 1;
            }
            c.swap(cn);
            if (store) tmp_levels.emplace_back(c.begin(), c.begin() + min(static_cast<int>(c.size()), n));
        }
        if (store) c_levels = move(tmp_levels);
        return p;
    }
    [[nodiscard]] int compare_suffix_pattern(int pos, const string &pat) const {
        int i = 0;
        while (pos + i < n and i < pat.size()) {
            if (s[pos + i] < pat[i]) return -1;
            if (s[pos + i] > pat[i]) return 1;
            i++;
        }
        if (i == pat.size()) return 0;
        return -1;
    }
    void build_kasai() {
        lcp.assign(max(0, n - 1), 0);
        int k = 0;
        for (int i = 0; i < n; i++) {
            const int pi = pos[i];
            if (pi == n - 1) { k = 0; continue; }
            const int j = suffix_array[pi + 1];
            while (i + k < n and j + k < n and s[i + k] == s[j + k]) ++k;
            lcp[pi] = k;
            if (k) --k;
        }
    }
    void build_rmq() {
        const int m = static_cast<int>(lcp.size());
        if (m == 0) { st.clear(); lg.assign(1, 0); return; }
        lg.assign(m + 1, 0);
        for (int i = 2; i <= m; ++i) lg[i] = lg[i >> 1] + 1;
        const int LOG = lg[m] + 1;
        st.assign(LOG, vector<int>(m));
        for (int i = 0; i < m; i++) st[0][i] = lcp[i];
        for (int k = 1; k < LOG; k++)
            for (int i = 0; i + (1 << k) <= m; ++i)
                st[k][i] = min(st[k - 1][i], st[k - 1][i + (1 << (k - 1))]);
    }
    [[nodiscard]] int rmq(const int l, const int r) const {
        if (l > r) return 0;
        const int len = r - l + 1;
        const int k = lg[len];
        return min(st[k][l], st[k][r - (1 << k) + 1]);
    }
};

struct SuffixAutomaton {
    struct state {
        ll cnt = 0;
        vector<int> inv_link;
        bool is_clone = false;
        unordered_map<char, int> next;
        int len = 0, link = -1, first_pos = -1;
    };

    int sz, last;
    vector<ll> dp;
    vector<state> st;
    bool kth_ready = false;
    bool occ_ready = false;
    vector<vector<pair<char,int>>> sorted_next;

    explicit SuffixAutomaton(const string &s) {
        st.reserve(s.size() << 1);
        st.emplace_back();
        st[0].len = 0, st[0].link = -1;
        sz = 1, last = 0;
        for (const char c : s) extend(c);
    }
    [[nodiscard]] ll count_different_substrings() const {
        ll res = 0;
        for (int i = 1; i < sz; i++) res += st[i].len - st[st[i].link].len;
        return res;
    }
    [[nodiscard]] ll total_length_diff_substings() const {
        ll res = 0;
        for(int i = 1; i < sz; i++) {
            const ll shortest = st[st[i].link].len + 1;
            const ll longest = st[i].len;
            const ll num_strings = longest - shortest + 1;
            const ll cur = num_strings * (longest + shortest) / 2;
            res += cur;
        }
        return res;
    }
    string kth_substring(ll k) {
        if (!kth_ready) prepare_kth();
        // if (!occ_ready) prepare_occurrence_info();
        if (k <= 0 or k > dp[0]) return {};
        int u = 0; string ans;
        while (k > 0)
        for (auto &[fst, snd] : sorted_next[u]) {
            const char c = fst;
            const int to = snd;
            if (const ll block = 1 + dp[to] /* ll block = dp[to] */ ; k > block) k -= block;
            else {
                ans.push_back(c);
                if (k == 1) return ans; // if (k <= st[to].cnt)
                k--, u = to; // k -= st[to].cnt
                break;
            }
        }
        return ans;
    }
    bool contains(const string &p) {
        int u = 0;
        for (char c : p) {
            if (!st[u].next.count(c)) return false;
            u = st[u].next[c];
        }
        return true;
    }
    ll occurrences_count(const string &p) {
        if (!occ_ready) prepare_occurrence_info();
        int u = 0;
        for (char c : p) {
            if (!st[u].next.count(c)) return 0;
            u = st[u].next.at(c);
        }
        return st[u].cnt;
    }
    int first_occurrence(const string &p) {
        if (!occ_ready) prepare_occurrence_info();
        int u = 0;
        for (char c : p) {
            if (!st[u].next.count(c)) return -1;
            u = st[u].next.at(c);
        }
        if (st[u].first_pos == -1) return -1;
        return st[u].first_pos - static_cast<int>(p.size()) + 1;
    }
    vector<int> all_occurrences(const string &p) {
        if (!occ_ready) prepare_occurrence_info();
        int u = 0;
        for (char c : p) {
            if (!st[u].next.count(c)) return {};
            u = st[u].next.at(c);
        }
        vector<int> res, stack;
        stack.reserve(64), stack.push_back(u);
        while (!stack.empty()) {
            const int v = stack.back(); stack.pop_back();
            if (!st[v].is_clone and st[v].first_pos != -1) {
                if (int start = st[v].first_pos - static_cast<int>(p.size()) + 1; start >= 0) res.push_back(start);
            }
            for (int to : st[v].inv_link) stack.push_back(to);
        }
        sort(res.begin(), res.end());
        res.erase(unique(res.begin(), res.end()), res.end());
        return res;
    }
    [[nodiscard]] string shortest_nonappearing(const string &alphabet) const {
        constexpr int INF = 1e9;
        int maxlen = 0;
        for (int i = 0; i < sz; i++) maxlen = max(maxlen, st[i].len);
        vector<int> cnt(maxlen + 1), order(sz);
        for (int i = 0; i < sz; i++) cnt[st[i].len]++;
        for (int i = 1; i <= maxlen; i++) cnt[i] += cnt[i - 1];
        for (int i = sz - 1; i >= 0; --i) order[--cnt[st[i].len]] = i;
        vector d(sz, INF);
        for (int idx = sz - 1; idx >= 0; --idx) {
            const int v = order[idx];
            bool missing = false;
            for (char c : alphabet) if (!st[v].next.count(c)) { d[v] = 1; missing = true; break; }
            if (missing) continue;
            int best = INF;
            for (char c : alphabet) best = min(best, d[st[v].next.at(c)]);
            if (best < INF) d[v] = 1 + best;
        }
        if (d[0] >= INF) return {};
        string ans;
        int v = 0;
        while (true) {
            bool done = false;
            for (char c : alphabet) if (!st[v].next.count(c)) { ans.push_back(c); done = true; break; }
            if (done) break;
            for (char c : alphabet) {
                if (const int to = st[v].next.at(c); d[to] == d[v] - 1) { ans.push_back(c); v = to; break; }
            }
        }
        return ans;
    }
    string longest_common_substring(const string &t) {
        int u = 0, l = 0, best = 0, bestpos = 0;
        for (int i = 0; i < t.size(); i++) {
            while (u and !st[u].next.count(t[i])) {
                u = st[u].link;
                l = st[u].len;
            }
            if (st[u].next.count(t[i])) {
                u = st[u].next[t[i]];
                l++;
            }
            if (l > best) {
                best = l;
                bestpos = i;
            }
        }
        return t.substr(bestpos - best + 1, best);
    }

private:
    void extend(const char c) {
        const int cur = sz++;
        int p = last;
        st.emplace_back();
        st[cur].len = st[last].len + 1;
        st[cur].first_pos = st[cur].len - 1;
        st[cur].cnt = 1;
        st[cur].is_clone = false;
        while (p != -1 and !st[p].next.count(c)) {
            st[p].next[c] = cur;
            p = st[p].link;
        }
        if (p == -1) st[cur].link = 0;
        else {
            if (const int q = st[p].next[c]; st[p].len + 1 == st[q].len) st[cur].link = q;
            else {
                const int clone = sz++;
                st.emplace_back();
                st[clone].len = st[p].len + 1;
                st[clone].next = st[q].next;
                st[clone].link = st[q].link;
                st[clone].is_clone = true;
                st[clone].first_pos = st[q].first_pos;
                st[clone].cnt = 0;
                while (p != -1 and st[p].next[c] == q) {
                    st[p].next[c] = clone;
                    p = st[p].link;
                }
                st[q].link = st[cur].link = clone;
            }
        }
        last = cur;
    }
    void prepare_kth() {
        if (kth_ready) return;
        sorted_next.assign(sz, {});
        for (int i = 0; i < sz; i++) {
            sorted_next[i].reserve(st[i].next.size());
            for (auto &kv : st[i].next) sorted_next[i].emplace_back(kv);
            sort(sorted_next[i].begin(), sorted_next[i].end(),
                 [](const pair<char,int> &a, const pair<char,int> &b){ return a.first < b.first; });
        }
        int maxlen = 0;
        for (int i = 0; i < sz; i++) maxlen = max(maxlen, st[i].len);
        vector<int> cnt(maxlen + 1), order(sz);
        for (int i = 0; i < sz; i++) cnt[st[i].len]++;
        for (int i = 1; i <= maxlen; i++) cnt[i] += cnt[i - 1];
        for (int i = sz - 1; i >= 0; --i) order[--cnt[st[i].len]] = i;
        dp.assign(sz, 0);
        for (int idx = sz - 1; idx >= 0; idx--) {
            const int u = order[idx];
            ll sum = 0; // ll sum = st[u].cnt
            for (const auto &val: st[u].next) sum += 1 + dp[val.second]; // sum += dp[v]
            dp[u] = sum;
        }
        kth_ready = true;
    }
    void prepare_occurrence_info() {
        if (occ_ready) return;
        int maxlen = 0;
        for (int i = 0; i < sz; i++) maxlen = max(maxlen, st[i].len);
        vector<int> cnt_len(maxlen + 1, 0), order(sz);
        for (int i = 0; i < sz; i++) ++cnt_len[st[i].len];
        for (int i = 1; i <= maxlen; i++) cnt_len[i] += cnt_len[i-1];
        for (int i = sz - 1; i >= 0; i--) order[--cnt_len[st[i].len]] = i;
        for (int idx = sz - 1; idx > 0; idx--) {
            const int v = order[idx];
            if (const int p = st[v].link; p != -1) st[p].cnt += st[v].cnt;
        }
        for (int v = 1; v < sz; v++)
            if (const int p = st[v].link; p != -1) st[p].inv_link.push_back(v);
        occ_ready = true;
    }
};
