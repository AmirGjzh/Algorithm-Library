#include <bits/stdc++.h>
using namespace std;
using ll = long long int;

/*============================================================================================================
Binomial Coefficients C(n,k) — multiple methods, ordered, CP usage guide

Description:
  • C(n,k) = n choose k = number of k-element subsets from n elements
  • Formula: C(n,k) = n! / (k! * (n-k)!)
  • Widely used in combinatorics, DP combinatorics, probability, counting problems

Implementations:

1) Multiplicative method with gcd reduction (mul_gcd)
   • Exact integer, reduces numerator/denominator by gcd before multiplication
   • Complexity: O(k), safe for moderate k and result fits in 64-bit
   • Use when n large, k small, want exact result

2) Improved incremental multiplicative method (mul_inc)
   • Incrementally multiplies while reducing by gcd, safer for larger k
   • Complexity: O(k²) worst-case due to gcd checks
   • Use when k moderate, need exact integer, want safer intermediate

3) Pascal triangle DP table (build_pascal)
   • Precompute full C[n][k] table for small max n
   • Complexity: O(n²) time and memory
   • Use when max n small (<= ~4000) and many queries needed

4) Factorials + inverse factorials modulo prime (FactMod)
   • Precompute fact[i] = i! % MOD, invfact[i] = (i!)^-1 % MOD
   • Complexity: O(MAXN) preprocess, O(1) query
   • Requires prime MOD
   • Use when many queries modulo a prime

5) Lucas theorem (Lucas)
   • Handles huge n (up to 1e18) with small prime modulus p
   • Decompose n, k in base p, multiply small-digit C(n_i,k_i) mod p
   • Complexity: O(log_p(n))
   • Use when n large, p small prime

6) Prime-power / CRT method (not implemented here)
   • Use for composite moduli or prime powers
   • Combine C(n,k) mod p^q via Legendre's formula and CRT
   • Complexity: problem-specific

Notes:
  • Overflow: multiplicative methods may exceed 64-bit for large k — careful
  • gcd reduction keeps intermediate small for exact computation
  • For modulo computations, Fermat's little theorem is used for inverses
  • For Lucas, need precomputed fact/ifact up to p-1
============================================================================================================*/

ll mul_gcd(const ll n, ll k) {
    if (k < 0 or k > n) return 0;
    k = min(k, n - k);
    if (k == 0) return 1;
    vector<ll> num(k), den(k);
    for (ll i = 0; i < k; i++) {
        num[i] = n - k + 1 + i;
        den[i] = i + 1;
    }
    for (ll i = 0; i < k; i++) {
        if (den[i] == 1) continue;
        for (ll j = 0; j < k and den[i] > 1; j++) {
            if (const ll g = gcd(den[i], num[j]); g > 1) {
                den[i] /= g;
                num[j] /= g;
            }
        }
    }
    ll res = 1;
    for (ll i = 0; i < k; ++i) res *= num[i];
    return res;
}

ll mul_inc(const ll n, ll k) {
    if (k < 0 or k > n) return 0;
    k = min(k, n - k);
    if (k == 0) return 1;
    vector<ll> den(k + 1);
    for (ll i = 1; i <= k; ++i) den[i] = i;
    ll res = 1;
    for (ll i = 1; i <= k; ++i) {
        ll cur = n - k + i;
        for (ll j = 1; j <= k and cur > 1; ++j) {
            if (const ll g = gcd(cur, den[j]); g > 1) {
                cur /= g;
                den[j] /= g;
            }
        }
        res *= cur;
        for (ll j = 1; j <= k; ++j) {
            if (den[j] == 1) continue;
            gcd(res % den[j], den[j]);
        }
    }
    return res;
}

vector<vector<ll>> build_pascal(const int N) {
    vector C(N + 1, vector<ll>(N + 1, 0));
    for (int n = 0; n <= N; ++n) {
        C[n][0] = C[n][n] = 1;
        for (int k = 1; k < n; ++k) C[n][k] = C[n - 1][k - 1] + C[n - 1][k];
    }
    return C;
}

struct FactMod {
    ll MOD;
    vector<ll> fact, invfact;

    explicit FactMod(const int N = 0, const ll mod = 1e9 + 7) {
        MOD = mod;
        if (N > 0) init(N);
    }
    void init(const int N) {
        fact.assign(N + 1, 1);
        invfact.assign(N + 1, 1);
        for (int i = 1; i <= N; ++i) fact[i] = fact[i - 1] * i % MOD;
        invfact[N] = mod_pow(fact[N], MOD - 2, MOD);
        for (int i = N; i >= 1; --i) invfact[i - 1] = invfact[i] * i % MOD;
    }
    static ll mod_pow(ll a, ll e, const ll mod) {
        ll r = 1 % mod;
        a %= mod;
        while (e) {
            if (e & 1) r = r * a % mod;
            a = a * a % mod;
            e >>= 1;
        }
        return r;
    }
    [[nodiscard]] ll C(const ll n, const ll k) const {
        if (k < 0 or k > n) return 0;
        return fact[n] * invfact[k] % MOD * invfact[n - k] % MOD;
    }
};

struct Lucas {
    static ll mod_prime_small(ll n, ll k, const ll p, const vector<ll> &fact, const vector<ll> &invfact) {
        if (k < 0 or k > n) return 0;
        ll res = 1;
        while (n > 0 or k > 0) {
            const ll ni = n % p;
            const ll ki = k % p;
            if (ki > ni) return 0;
            ll cur = fact[ni] * invfact[ki] % p;
            cur = cur * invfact[ni - ki] % p;
            res = res * cur % p;
            n /= p, k /= p;
        }
        return res;
    }
    static ll mod_pow(ll a, ll e, const ll mod) {
        ll r = 1 % mod;
        a %= mod;
        while (e) {
            if (e & 1) r = r * a % mod;
            a = a * a % mod;
            e >>= 1;
        }
        return r;
    }
    static pair<vector<ll>, vector<ll>> build_fact_inv_for_p(const int p) {
        vector<ll> fact(p), invfact(p);
        fact[0] = 1;
        for (int i = 1; i < p; ++i) fact[i] = fact[i - 1] * i % p;
        invfact[p - 1] = mod_pow(fact[p - 1], p - 2, p);
        for (int i = p - 1; i >= 1; --i) invfact[i - 1] = invfact[i] * i % p;
        return {fact, invfact};
    }
};
