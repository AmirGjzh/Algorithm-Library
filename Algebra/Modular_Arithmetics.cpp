#include <bits/stdc++.h>
using namespace std;
using ll = long long int;

/*============================================================================================================
Modular Arithmetic & Number Theory Toolkit

Overview:
  This file provides a collection of classic number-theoretic algorithms
  used in modular arithmetic, cryptography, and competitive programming.
  It covers modular inverses, the Chinese Remainder Theorem, discrete logarithms,
  primitive roots, and discrete k-th roots modulo a prime.

1. Modular Inverse:
  - Finds x such that:
      a·x ≡ 1 (mod m)
  - Exists iff gcd(a, m) = 1
  - Implemented using:
      • Extended Euclidean Algorithm (general m)
      • Recursive formula for prime modulus
  - Also supports computing inverses for all values 1 … m−1 in O(m) time
      when m is prime.

2. Chinese Remainder Theorem (CRT):
  - Solves systems of congruences:
      x ≡ a_i (mod m_i)
  - Requires moduli m_i to be pairwise coprime
  - If moduli are not coprime, congruences must be merged carefully
      and consistency checked via gcd conditions
  - Returns the unique solution modulo M = ∏ m_i

3. Discrete Logarithm (Baby-Step Giant-Step):
  - Finds the smallest x ≥ 0 such that:
      a^x ≡ b (mod m)
  - Handles cases where gcd(a, m) ≠ 1 by reducing the problem
  - Uses a hash map for baby steps and modular exponentiation for giant steps

4. Primitive Roots:
  - A number g is a primitive root modulo p if its powers generate
      all nonzero residues module p
  - For a prime p, primitive roots always exist
  - A value g is a primitive root iff:
      g^(φ(p) / q) ≠ 1 (mod p) for all prime divisors q of φ(p)
  - This implementation finds the smallest primitive root for prime p

5. Discrete k-th Roots Modulo a Prime:
  - Solves:
      x^k ≡ a (mod n), where n is prime
  - Uses primitive roots and a Baby-Step Giant-Step-style reduction
  - Returns all valid solutions in the range [0 … n−1]

Conventions:
  • Moduli are assumed positive
  • All modular inverses return −1 if no inverse exists
  • Input vectors for CRT use pairs {remainder, modulus}
  • All solutions are normalized modulo the product or modulus

Complexity Summary:
  • Modular inverse (single):           O(log(m))
  • Inverses 1 … m−1 (prime m):         O(m)
  • Chinese Remainder Theorem:          O(k + log(M))
  • Discrete logarithm (BSGS):          O(√m·log(m))
  • Primitive root (prime p):           O(factors(φ(p))·log(p))
  • Discrete k-th root (prime n):       O(√n·log(n))

Applications:
  • Cryptography (RSA, Diffie–Hellman, discrete log problems)
  • Modular equation solving
  • Competitive programming
  • Number-theoretic algorithms and proofs

Notes:
  • CRT requires pairwise coprime moduli for correctness
  • Discrete logarithms may not exist for all inputs
  • Floating-point operations are not usedm all math is exact
  • Primitive roots are only searched for prime moduli here
============================================================================================================*/

ll mod_inverse(ll a, const ll m) {
    a %= m;
    if (a < 0) a += m;
    ll b = m, x0 = 1, x1 = 0;
    while (b) {
        const ll q = a / b;
        tie(a, b) = make_pair(b, a - q * b);
        tie(x0, x1) = make_pair(x1, x0 - q * x1);
    }
    if (a != 1) return -1;
    if (x0 < 0) x0 += m;
    return x0;
}

ll modular_inverse(ll a, const ll m) {
    a %= m;
    return a <= 1 ? a : m - m / a * modular_inverse(m % a, m) % m;
}

vector<ll> modular_inverse_1_m(const ll m) {
    vector<ll> inv(m);
    inv[1] = 1;
    for (int i = 2; i < m; i++) inv[i] = m - m / i * inv[m % i] % m;
    return inv;
}

ll chinese_remainder_theorem(const vector<pair<ll, ll>> &congruences) {
    ll M = 1, solution = 0;
    for (const auto &val : congruences) M *= val.second;
    for (const auto &[fst, snd] : congruences) {
        const ll a_i = fst;
        const ll M_i = M / snd;
        const ll N_i = modular_inverse(M_i, snd);
        solution = (solution + a_i * M_i % M * N_i) % M;
    }
    return solution;
}

ll baby_step_giant_step(ll a, ll b, ll m) {
    a %= m, b %= m;
    if (a == 0) return b == 0 ? 1 : -1;
    ll k = 1, add = 0, g;
    while ((g = __gcd(a, m)) > 1) {
        if (b == k) return add;
        if (b % g) return -1;
        b /= g, m /= g, add++, k = (k * a / g) % m;
    }
    const ll n = static_cast<ll>(sqrt(m)) + 1;
    ll an = 1;
    for (ll i = 0; i < n; i++) an = an * a % m;
    unordered_map<ll, ll> vals;
    for (ll q = 0, cur = b; q <= n; q++) vals[cur] = q, cur = cur * a % m;
    for (ll p = 1, cur = k; p <= n; p++) {
        cur = cur * an % m;
        if (vals.count(cur)) return n * p - vals[cur] + add;
    }
    return -1;
}

ll binpow(ll a, ll b, const ll mod) {
    a %= mod;
    ll res = 1;
    while (b > 0) { if (b & 1) res = res * a % mod; a = a * a % mod, b >>= 1; }
    return res;
}

ll Phi(ll n) {
    ll result = n;
    for (ll i = 2; i * i <= n; i++) if (n % i == 0) {
        while (n % i == 0) n /= i;
        result -= result / i;
    }
    if (n > 1) result -= result / n;
    return result;
}

ll find_primitive_root(const ll p) {
    vector<ll> fact;
    const ll phi = Phi(p);
    ll n = phi;
    for (ll i = 2; i * i <= n; i++) if (n % i == 0) {
        fact.push_back(i);
        while (n % i == 0) n /= i;
    }
    if (n > 1) fact.push_back(n);
    for (ll res = 2; res <= p; res++) {
        if (__gcd(res, p) != 1) continue;
        bool ok = true;
        for (size_t i = 0; i < fact.size() && ok; i++)
            ok &= binpow(res, phi / fact[i], p) != 1;
        if (ok) return res;
    }
    return -1;
}

vector<ll> find_all_discrete_roots(const ll n, const ll k, const ll a) {
    if (a == 0) return {0};
    const ll g = find_primitive_root(n);
    const ll sq = static_cast<ll>(sqrt(n)) + 1;
    vector<pair<ll, ll>> dec(sq);
    for (ll i = 1; i <= sq; i++)
        dec[i - 1] = {binpow(g, i * sq * k % (n - 1), n), i};
    sort(dec.begin(), dec.end());
    ll any_ans = -1;
    for (ll i = 0; i < sq; i++) {
        ll my = binpow(g, i * k % (n - 1), n) * a % n;
        if (auto it = lower_bound(dec.begin(), dec.end(), make_pair(my, 0LL)); it != dec.end() && it->first == my) {
            any_ans = it->second * sq - i;
            break;
        }
    }
    if (any_ans == -1) return {};
    const ll delta = (n - 1) / __gcd(k, n - 1);
    vector<ll> ans;
    for (ll cur = any_ans % delta; cur < n - 1; cur += delta)
        ans.push_back(binpow(g, cur, n));
    sort(ans.begin(), ans.end());
    return ans;
}
