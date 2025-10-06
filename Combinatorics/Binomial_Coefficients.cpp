#include <bits/stdc++.h>
using namespace std;
using ll = long long int;
using i128 = __int128_t;

/*============================================================================================================
Binomial Coefficients - summary + multiple implementations

Description:
  • Binomial coefficient C(n, k) = n choose k = number of k-element subsets from n elements
  • Analytic formula: C(n,k) = n!/(k!.(n-k)!)
  • Many algorithmic variants are useful depending on constraints (n, k sizes, modulo, many queries)
  • Common competitive-programming implementations:
  
    1) Multiplicative integer method (exact integer, reduced by gcds) - good for single nCr without overflow:
        Complexity: O(k) where k = min(k, n-k)
        Strategy:   keep arrays of numerators and denominators, reduce by gcd to avoid overflow,
                    multiply using 128-bit intermediate to be safer
        Returns:    exact integer C(n,k). Works up to the limits where result fits in signed 64-bit (or more if you change return type)
    2) Improved multiplicative using __int128 for safer multiplication:
        This version multiplies numerator/denominator incrementally while reducing by gcd at each step
        Safe and fast for moderate k, uses __int128 for intermediate products
    3) Pascal triangle DP (precompute full table) - O(N^2) preprocessing, O(1) queries:
        Use when max n is small (e.g., n <= 2000..5000 depending on memory/time)
        Complexity: O(N^2) time and memory to build whole triangle
    4) Precompute factorials + inverse factorials modulo prime (O(1) queries after O(N) preprocess):
        Precompute fact[i] = i! % MOD and invfact[i] = (i!)^-1 % MOD for i up to MAXN
        Then mod(n,k) = fact[n] * invfact[k] % MOD * invfact[n-k] % MOD
        Requirements: MOD must be prime (so inverse via Fermat's little theorem)
        Precompute complexity: O(MAXN) time, each query O(1)
    5) Lucas theorem for large n with small prime modulus:
        When n is large (e.g., up to 1e18) but modulus p is small prime (e.g., p <= 1e6),
        Lucas reduces to digits in base p:
            C(n,k) mod p = Prod_{i} C(n_i, k_i) mod p
        where n_i, k_i are digits of n and k in base p
        Requires precomputed small factorials mod p for 0..p-1
    6) Notes & references for modulo prime-power and arbitrary modulus (more complex; see referenced page):
        If modulus is a prime power p^q, computing C(n,k) mod p^q is more involved:
          - you must track multiplicity of prime p in n!, k!, (n-k)! (Legendre's formula),
            compute factorials with p factors removed modulo p^q, and combine using multiplicative inverses
            modulo p^q for non-p factors
        For arbitrary modulus m, decompose m into prime powers, compute C(n,k) mod p^q for each component,
        then combine results via Chinese Remainder Theorem (CRT)

  • Use each method depending on limits:
    - Small n (<= few thousands): Pascal table or factorial array
    - Large n but small k: multiplicative method (O(k)) is best
    - Many queries with modulo prime: precompute factorials + invfacts
    - Large n and small prime modulus: Lucas
    - Modulo ≠ prime or prime power: see CP-algorithms section for prime-power decomposition + CRT
Applications:
  • Combinatorics, probability, counting problems, polynomial expansions, DP combinatorics

Notes:
  • For modulo p (prime) use Fermat inverse (a^(p-2) mod p)
  • For composite modulo you generally need to decompose into prime powers and use Chinese Remainder Theorem
============================================================================================================*/

ll mul_gcd(ll n, ll k) {
    if (k < 0 || k > n) return 0;
    k = min(k, n - k);
    if (k == 0) return 1;
    vector<ll> num(k), den(k);
    for (ll i = 0; i < k; i++) {
        num[i] = n - k + 1 + i; 
        den[i] = i + 1;        
    }
    for (ll i = 0; i < k; i++) {
        if (den[i] == 1) continue;
        for (ll j = 0; j < k && den[i] > 1; j++) {
            ll g = gcd(den[i], num[j]);
            if (g > 1) {
                den[i] /= g;
                num[j] /= g;
            }
        }
    }
    i128 res = 1;
    for (ll i = 0; i < k; ++i) res *= (i128)num[i];
    return (ll) res;
}

ll mul_inc(ll n, ll k) {
    if (k < 0 || k > n) return 0;
    k = min(k, n - k);
    if (k == 0) return 1;
    vector<ll> den(k + 1);
    for (ll i = 1; i <= k; ++i) den[i] = i;
    i128 res = 1;
    for (ll i = 1; i <= k; ++i) {
        i128 cur = n - k + i;
        for (ll j = 1; j <= k && cur > 1; ++j) {
            ll g = gcd((ll) cur, den[j]);
            if (g > 1) {
                cur /= g;
                den[j] /= g;
            }
        }
        res *= cur;
        for (ll j = 1; j <= k; ++j) {
            if (den[j] == 1) continue;
            ll g = (ll) gcd((ll) (res % den[j]), den[j]);
        }
    }
    return (ll) res;
}

vector<vector<ll>> build_pascal(int N) {
    vector<vector<ll>> C(N + 1, vector<ll>(N + 1, 0));
    for (int n = 0; n <= N; ++n) {
        C[n][0] = C[n][n] = 1;
        for (int k = 1; k < n; ++k) C[n][k] = C[n - 1][k - 1] + C[n - 1][k];
    }
    return C;
}

struct FactMod {
    ll MOD;
    vector<ll> fact, invfact;
    FactMod(int N = 0, ll mod = 1'000'000'007) {
        MOD = mod;
        if (N > 0) init(N);
    }
    void init(int N) {
        fact.assign(N + 1, 1);
        invfact.assign(N + 1, 1);
        for (int i = 1; i <= N; ++i) fact[i] = ((__int128)fact[i - 1] * i) % MOD;
        invfact[N] = mod_pow(fact[N], MOD - 2, MOD);
        for (int i = N; i >= 1; --i) invfact[i - 1] = ((__int128)invfact[i] * i) % MOD;
    }
    ll mod_pow(ll a, ll e, ll mod) {
        ll r = 1 % mod;
        a %= mod;
        while (e) {
            if (e & 1) r = ((__int128) r * a) % mod;
            a = ( (__int128) a * a) % mod;
            e >>= 1;
        }
        return r;
    }
    ll C(ll n, ll k) {
        if (k < 0 || k > n) return 0;
        return ((__int128) fact[n] * invfact[k] % MOD) * invfact[n - k] % MOD;
    }
};

struct Lucas {
    ll mod_prime_small(ll n, ll k, ll p, const vector<ll>& fact, const vector<ll>& invfact) {
        if (k < 0 || k > n) return 0;
        ll res = 1;
        while (n > 0 || k > 0) {
            ll ni = n % p;
            ll ki = k % p;
            if (ki > ni) return 0;
            ll cur = ((__int128)fact[ni] * invfact[ki]) % p;
            cur = ((__int128)cur * invfact[ni - ki]) % p;
            res = ((__int128)res * cur) % p;
            n /= p;
            k /= p;
        }
        return res;
    }
    ll mod_pow(ll a, ll e, ll mod) {
        ll r = 1 % mod;
        a %= mod;
        while (e) {
            if (e & 1) r = ((__int128) r * a) % mod;
            a = ( (__int128) a * a) % mod;
            e >>= 1;
        }
        return r;
    }
    pair<vector<ll>, vector<ll>> build_fact_inv_for_p(int p) {
        vector<ll> fact(p), invfact(p);
        fact[0] = 1;
        for (int i = 1; i < p; ++i) fact[i] = ( (__int128)fact[i - 1] * i ) % p;
        invfact[p - 1] = mod_pow(fact[p - 1], p - 2, p);
        for (int i = p - 1; i >= 1; --i) invfact[i - 1] = ((__int128)invfact[i] * i) % p;
        return {fact, invfact};
    }
};
