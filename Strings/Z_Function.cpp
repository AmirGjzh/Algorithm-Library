#include <bits/stdc++.h>
using namespace std;
using ll = long long int;

/*============================================================================================================
Z-Function Algorithm and Applications

1. Z-Function (z[i])
   - Length of longest substring starting at i that matches prefix of s
   - Maintains [l, r] interval of rightmost matched segment
   - Each character compared at most twice → O(n) complexity

2. Pattern Matching with Z-Function
   - Similar to KMP but uses Z-values for prefix matching
   - Pattern + '#' + text construction for separate matching

3. Applications of Z-Function
   - Counting distinct substrings in O(n²)
   - Finding shortest repeating substring
   - Finding palindromic prefixes
   - Counting prefix occurrences in string
   - Wildcard pattern matching

Key Insights:
  • When i ≤ r, reuse: z[i] = min(r-i+1, z[i-l])
  • Update [l, r] when i+z[i]-1 > r
  • Similar to prefix function but measures matches from start position

Complexity Summary:
  • Z-function computation: O(n) time, O(n) space
  • Pattern matching: O(n+m) time, O(n+m) space
  • Distinct substrings: O(n²) time, O(n) space

Conventions:
  • All indices are 0-based, z[0] = 0 (or n, but 0 is standard)
  • Separator '#' must not appear in pattern or text
  • Palindromic prefixes: check z[n + 1 + (n - i)] ≥ i for prefix length i
  • Wildcard matching: '?' matches any character

Comparison with Prefix Function:
  • Z-function: prefix match lengths starting at each position
  • Prefix function: border lengths ending at each position
  • Both solve similar problems with different formulations
============================================================================================================*/

vector<int> z_function(const string &s) {
    int l = 0, r = 0;
    const int n = static_cast<int>(s.size());
    vector z(n, 0);
    for (int i = 1; i < n; i++) {
        if (i <= r) z[i] = min(r - i + 1, z[i - l]);
        while (i + z[i] < n and s[z[i]] == s[i + z[i]]) z[i]++;
        if (i + z[i] - 1 > r) l = i, r = i + z[i] - 1;
    }
    return z;
}

vector<int> z_search(const string &p, const string &t) {
    if (p.empty()) return {};
    const string s = p + '#' + t;
    const int n = static_cast<int>(s.size()), m = static_cast<int>(t.size());
    const vector<int> z = z_function(s);
    vector<int> matches;
    for (int i = n + 1; i < n + m + 1; i++) if (z[i] == n) matches.push_back(i - n - 1);
    return matches;
}

int count_distinct_substrings_z(const string &s) {
    const int n = static_cast<int>(s.size());
    int total = 0;
    string cur;
    for (int i = 0; i < n; i++) {
        cur += s[i];
        string rev = cur;
        reverse(rev.begin(), rev.end());
        vector<int> z = z_function(rev);
        const int z_max = *max_element(z.begin(), z.end());
        total += static_cast<int>(cur.size()) - z_max;
    }
    return total;
}

string shortest_repeating_substring_z(const string &s) {
    const int n = static_cast<int>(s.size());
    const vector<int> z = z_function(s);
    for (int i = 1; i < n; i++) if (n % i == 0 and i + z[i] == n) return s.substr(0, i);
    return s;
}

vector<int> find_palindromic_prefixes(const string &s) {
    const int n = static_cast<int>(s.size());
    string rev = s;
    reverse(rev.begin(), rev.end());
    const string t = s + '#' + rev;
    const vector<int> z = z_function(t);
    vector<int> palindromes;
    for (int i = 1; i <= n; i++) if (z[n + 1 + (n - i)] >= i) palindromes.push_back(i);
    return palindromes;
}

vector<int> count_prefix_occurrences_z(const string &s) {
    const int n = static_cast<int>(s.size());
    const vector<int> z = z_function(s);
    vector diff(n + 2, 0);
    for (int i = 1; i < n; i++) if (z[i] > 0) diff[1]++, diff[z[i] + 1]--;
    vector ans(n + 1, 0);
    int current = 0;
    for (int i = 1; i <= n; i++) {
        current += diff[i];
        ans[i] = current;
    }
    for (int i = 1; i <= n; i++) ans[i]++;
    return ans;
}

vector<int> z_search_with_wildcard(const string &p, const string &t) {
    const int n = static_cast<int>(p.size()), m = static_cast<int>(t.size());
    auto char_match = [](const char a, const char b) -> bool {
        if (a == '#' or b == '#') return false;
        return a == '?' || b == '?' || a == b;
    };
    auto wildcard_z_function = [&](const string &s) -> vector<int> {
        int l = 0, r = 0;
        const int len = static_cast<int>(s.size());
        vector z(len, 0);
        for (int i = 1; i < len; i++) {
            if (i <= r) z[i] = min(r - i + 1, z[i - l]);
            while (i + z[i] < len && char_match(s[z[i]], s[i + z[i]])) z[i]++;
            if (i + z[i] - 1 > r) l = i, r = i + z[i] - 1;
        }
        return z;
    };
    const string s = p + '#' + t;
    const vector<int> z = wildcard_z_function(s);
    vector<int> matches;
    for (int i = n + 1; i < n + 1 + m; i++) if (z[i] == n) matches.push_back(i - n - 1);
    return matches;
}
