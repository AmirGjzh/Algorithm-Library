#include <bits/stdc++.h>
using namespace std;
using ll = long long int;

/*============================================================================================================
Hashing & Pattern Matching

1. String Hashing (Polynomial Rolling Hash)
   - Reverse hash formula: hash(s[0..i]) = Σ s[k] * p^(i-k) mod m
   - Enables O(1) substring queries after O(n) preprocessing
   - Double hashing reduces collision probability to ~10⁻¹⁸
   - Applications: pattern matching, substring counting, palindrome checking

2. Double Hashing
   - Uses two independent hash functions for collision resistance
   - Returns tuple: (hash1, hash2, length) for complete substring ID

3. Rabin-Karp Pattern Matching
   - O(n+m) average pattern matching using rolling hashes
   - Sliding window with O(1) hash comparison
   - Requires hash verification for absolute certainty

Conventions:
  • All indices are 0-based
  • Substring [l, r] includes both endpoints
  • Empty patterns are handled as edge cases
  • Uses (unsigned char) for extended ASCII support

Complexity Summary:
  • StringHasher construction: O(n) time, O(n) space
  • get_hash query: O(1) time
  • DoubleHasher: Same as StringHasher ×2
  • Rabin-Karp: O(n+m) time, O(n+m) space

Common Parameters:
  • p (base): 31 (lowercase), 53 (case-sensitive), 127 (ASCII)
  • m (modulus): 1e9+7, 1e9+9 (primes) or 2^64 via unsigned overflow
  • Default: p=31, m=1e9+9 for single hashing

Critical Notes:
  • Character mapping uses (unsigned char) for 0-255 range
  • Power array must be precomputed up to maximum substring length
  • Hash collisions are probabilistic - use DoubleHasher for critical apps
============================================================================================================*/

struct StringHasher {
    int p, m;
    vector<ll> pref_hash, p_pow;

    explicit StringHasher(const string &s, const int p = 31, const int m = 1e9 + 9) : p(p), m(m) {
        const int n = static_cast<int>(s.size());
        p_pow.assign(n + 1, 1);
        pref_hash.assign(n + 1, 0);
        for (int i = 1; i <= n; i++) p_pow[i] = p_pow[i - 1] * 1LL * p % m;
        for (int i = 0; i < n; i++) pref_hash[i + 1] = (pref_hash[i] * p + static_cast<unsigned char>(s[i])) % m;
    }
    [[nodiscard]] ll get_hash(const int l, const int r) const {
        ll hash_val = (pref_hash[r + 1] - pref_hash[l] * p_pow[r - l + 1]) % m;
        if (hash_val < 0) hash_val += m;
        return hash_val;
    }
};

struct DoubleHasher {
    StringHasher h1, h2;

    explicit DoubleHasher(const string &s) :
        h1(s, 131, 1e9 + 7),
        h2(s, 257, 1e9 + 9) {}
    [[nodiscard]] tuple<ll, ll, int> get_hash(const int l, const int r) const {
        return {h1.get_hash(l, r), h2.get_hash(l, r), r - l + 1};
    }
};

ll calc_hash(const string &s) {
    constexpr int p = 31, m = 1e9 + 9;
    ll hash_value = 0;
    for (const char c : s) hash_value = (hash_value * p + static_cast<unsigned char>(c)) % m;
    return hash_value;
}

vector<int> rabin_karp(const string &s, const string &t) {
	const int S = static_cast<int>(s.size()), T = static_cast<int>(t.size());
	const StringHasher sh(s), th(t);
	vector<int> pos;
    const ll pattern_hash = sh.get_hash(0, S - 1);
	for (int i = 0; i + S - 1 < T; i++)
		if (pattern_hash == th.get_hash(i, i + S - 1)) pos.push_back(i);
	return pos;
}
















/* ==========================
Number hashing for probabilty problems
like even occurnece check of all numbers in
an interval with XOR calculation
========================== */

uint64_t splitmix64(uint64_t x) {
    x += 0x9e3779b97f4a7c15;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9;
    x = (x ^ (x >> 27)) * 0x94d049bb133111eb;
    return x ^ (x >> 31);
}
