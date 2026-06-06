#include <bits/stdc++.h>
using namespace std;
using ll = long long int;

/*============================================================================================================
Prefix Function & Knuth-Morris-Pratt Algorithm

1. Prefix Function (π[i])
   - Length of longest proper prefix of s[0..i] that is also a suffix
   - Computed in O(n) time using iterative matching with fallback

2. Knuth-Morris-Pratt (KMP) Pattern Matching
   - O(n+m) pattern matching by avoiding redundant comparisons
   - Uses prefix function to skip unnecessary character checks

3. Applications of Prefix Function
   - Counting prefix occurrences in same string or text
   - Counting distinct substrings in O(n²)
   - Finding shortest repeating substring
   - Building prefix automaton for streaming pattern matching

Key Concepts:
  • π[i] ≤ π[i-1] + 1 (prefix function increases by at most 1)
  • When mismatch occurs, fallback to next best candidate using π
  • Pattern + '#' + text concatenation separates pattern from text

Complexity Summary:
  • Prefix function: O(n) time, O(n) space
  • KMP search: O(n+m) time, O(n+m) space
  • Prefix automaton: O(26n) time, O(26n) space
  • Distinct substrings: O(n²) time, O(n) space

Conventions:
  • All indices are 0-based
  • π[0] = 0 (empty string)
  • Separator '#' must not appear in pattern or text
  • Automaton assumes lowercase alphabet (a-z)

Important Notes:
  • Proper prefix excludes the whole string itself
  • For pattern matching, check when π[i] == pattern length
  • When counting prefixes, each prefix counts itself once
  • Shortest repeating substring length = n - π[n-1]
============================================================================================================*/

vector<int> prefix_function(const string &s) {
    const int n = static_cast<int>(s.size());
    vector pi(n, 0);
    for (int i = 1; i < n; i++) {
        int j = pi[i - 1];
        while (j > 0 and s[i] != s[j]) j = pi[j - 1];
        if (s[i] == s[j]) j++;
        pi[i] = j;
    }
    return pi;
}

vector<int> KMP(const string &s, const string &t) {
    const int n = static_cast<int>(s.size()), m = static_cast<int>(t.size());
    const string st = s + '#' + t;
    const vector<int> pi = prefix_function(st);
    vector<int> pos;
    for (int i = n + 1; i < n + m + 1; i++) if (pi[i] == n) pos.push_back(i - 2 * n);
    return pos;
}

vector<int> count_prefix_occurrences(const string &s) {
    const int n = static_cast<int>(s.size());
    vector ans(n + 1, 0);
    const vector<int> pi = prefix_function(s);
    for (int i = 0; i < n; i++) ans[pi[i]]++;
    for (int i = n - 1; i > 0; i--) ans[pi[i - 1]] += ans[i];
    for (int i = 0; i <= n; i++) ans[i]++;
    return ans;
}

vector<int> count_prefix_occurrences_in_text(const string &s, const string &t) {
    const int n = static_cast<int>(s.size()), m = static_cast<int>(t.size());
    const string st = s + '#' + t;
    vector ans(n + 1, 0);
    const vector<int> pi = prefix_function(st);
    for (int i = n + 1; i < n + 1 + m; i++) ans[pi[i]]++;
    for (int i = n - 1; i > 0; i--) ans[pi[i - 1]] += ans[i];
    for (int i = 0; i <= n; i++) ans[i]++;
    return ans;
}

int count_distinct_substrings(const string &s) {
    const int n = static_cast<int>(s.size());
    int total = 0;
    string cur;
    for (int i = 0; i < n; i++) {
        cur += s[i];
        string rev = cur;
        reverse(rev.begin(), rev.end());
        vector<int> pi = prefix_function(rev);
        const int pi_max = *max_element(pi.begin(), pi.end());
        total += static_cast<int>(cur.size()) - pi_max;
    }
    return total;
}

string shortest_repeating_substring(const string &s) {
    const int n = static_cast<int>(s.size());
    const vector<int> pi = prefix_function(s);
    if (const int len = n - pi[n - 1]; n % len == 0) return s.substr(0, len);
    return s;
}

vector<vector<int>> build_prefix_automaton(string s) {
    s += '#';
    const int n = static_cast<int>(s.size());
    const vector<int> pi = prefix_function(s);
    vector aut(n, vector(26, 0));
    for (int i = 0; i < n; i++)
        for (int c = 0; c < 26; c++)
            if (i > 0 and 'a' + c != s[i]) aut[i][c] = aut[pi[i - 1]][c];
            else aut[i][c] = i + ('a' + c == s[i]);
    return aut;
}
