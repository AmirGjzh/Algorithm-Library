#include <bits/stdc++.h>
using namespace std;
using u64 = uint64_t;
using u128 = __uint128_t;
using ll = long long int;
constexpr int MOD = 1e9 + 7;
static mt19937_64 rng(random_device{}());

/*============================================================================================================
Primality Testing & Integer Factorization Utilities

Overview:
  This file contains functions for primality checking, sieves, and integer factorization,
  covering both deterministic and probabilistic approaches. It also includes
  advanced factorization using Pollard’s Rho combined with Miller–Rabin.

1. Basic Primality Test (Trial Division)
  - Check divisibility up to √n
  - Complexity: O(√n)
  - Usage: small integers, simple checks

2. Probabilistic Primality Tests
  - Fermat Test: based on a^(n-1) ≡ 1 (mod n), k iterations
  - Miller–Rabin Test:
      • Decompose n-1 = d·2^s
      • Test random bases or fixed deterministic bases for 64-bit integers
  - Complexity: O(k·log³(n)) per test
  - Usage: large integers where deterministic trial division is slow

3. Sieve of Eratosthenes
  - Generates primes ≤ n in O(n·log(log(n)))
  - Returns boolean array is_prime[0…n]

4. Segmented Sieve
  - Generates primes in range [L, R] using primes ≤ √R
  - Complexity: O((R-L)·log(log(R)) + √R)

5. Linear Sieve
  - Computes lp[i] = least prime divisor for all i ≤ n in total O(n)
  - Also generates all primes up to n

6. Trial Division Factorization
  - Extracts prime factors of n in O(√n)
  - Appends remaining n if > 1

7. Advanced Factorization (Pollard’s Rho + Miller–Rabin)
  - Miller–Rabin: deterministic 64-bit bases for primality check
  - Pollard’s Rho: probabilistic factorization of composites
  - Recursively decomposes n into primes
  - Expected complexity per factor: O(√√n)
  - Applications:
      • Cryptography (weak RSA key analysis)
      • Competitive programming requiring full prime lists
      • Large integer factorization in number-theory problems

Conventions:
  • Random number generation via mt19937_64
  • Deterministic bases chosen for Miller–Rabin 64-bit reliability
============================================================================================================*/

bool is_prime(const ll n) {
    for (ll i = 2; i * i <= n; i++) if (n % i == 0) return false;
    return n >= 2;
}

struct ProbablyPrime {
    static u64 binpow(u64 a, u64 b, const u64 mod) {
        a %= mod;
        u64 res = 1;
        while (b > 0) {
            if (b & 1) res = static_cast<u128>(res) * a % mod; a = static_cast<u128>(a) * a % mod; b >>= 1;
        }
        return res;
    }
    static bool probably_prime_fermat(const int n, const int iter = 5) {
        if (n < 4) return n == 2 or n == 3;
        for (int i = 0; i < iter; i++) {
            uniform_int_distribution<u64> dist(2, n - 2);
            if (const u64 a = dist(rng); binpow(a, n - 1, n) != 1) return false;
        }
        return true;
    }
    static bool check_composite(const u64 n, const u64 a, const u64 d, const int s) {
        u64 x = binpow(a, d, n);
        if (x == 1 or x == n - 1) return false;
        for (int r = 1; r < s; r++) { x = static_cast<u128>(x) * x % n; if (x == n - 1) return false; }
        return true;
    }
    static bool miller_rabin1(const u64 n, const int iter = 5) {
        if (n < 4) return n == 2 or n == 3;
        u64 s = 0, d = n - 1;
        while ((d & 1) == 0) d >>= 1, s++;
        for (int i = 0; i < iter; i++) {
            uniform_int_distribution<u64> dist(2, n - 2);
            if (const u64 a = dist(rng); check_composite(n, a, d, static_cast<int>(s))) return false;
        }
        return true;
    }
    static bool miller_rabin2(const u64 n) {
        if (n < 2) return false;
        u64 r = 0, d = n - 1;
        while ((d & 1) == 0) d >>= 1, r++;
        for (const int a : {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37}) {
            if (n == a) return true;
            if (check_composite(n, a, d, static_cast<int>(r))) return false;
        }
        return true;
    }
};

vector<bool> sieve_of_eratosthenes(const int n) {
    vector is_prime(n + 1, true);
    is_prime[0] = is_prime[1] = false;
    for (int i = 2; i * i <= n; i++)
        if (is_prime[i]) for (int j = i * i; j <= n; j += i) is_prime[j] = false;
    return is_prime;
}

vector<char> segmented_sieve(const int L, const int R) {
    const int lim = static_cast<int>(sqrt(R));
    vector<int> primes;
    vector<char> mark(lim + 1, false);
    for (int i = 2; i <= lim; i++)
        if (!mark[i]) {
            primes.emplace_back(i);
            for (int j = i * i; j <= lim; j += i) mark[j] = true;
        }
    vector<char> is_prime(R - L + 1, true);
    for (const int i : primes)
        for (int j = max(i * i, (L + i - 1) / i * i); j <= R; j += i) is_prime[j - L] = false;
    if (L == 1) is_prime[0] = false;
    return is_prime;
}

vector<int> linear_sieve(const int n) {
    vector<int> primes;
    vector<int> lp(n + 1);
    for (int i = 2; i <= n; i++) {
        if (lp[i] == 0) lp[i] = i, primes.push_back(i);
        for (int j = 0; i * primes[j] <= n; j++) { lp[i * primes[j]] = primes[j]; if (primes[j] == lp[i]) break; }
    }
    return lp;
}

vector<ll> trial_division(ll n) {
    vector<ll> factorization;
    for (int d = 2; d * d <= n; d++) while (n % d == 0) {factorization.push_back(d); n /= d;}
    if (n > 1) factorization.push_back(n);
    return factorization;
}

struct Factorization {
    map<ll, int> factors;

    static ll mul_mod(const ll a, const ll b, const ll m) {
        return static_cast<ll>(static_cast<u128>(a) * b % m);
    }
    static ll binpow(ll a, ll b, const ll m) {
        ll res = 1;
        a %= m;
        while (b > 0) { if (b & 1) res = mul_mod(res, a, m); a = mul_mod(a, a, m); b >>= 1; }
        return res;
    }
    static bool is_prime(const ll n) {
        if (n < 2) return false;
        for (const ll p : { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37 }) if (n % p == 0) return n == p;
        ll d = n - 1, s = 0;
        while ((d & 1) == 0) d >>= 1, s++;
        for (const u64 a : { 2ULL, 325ULL, 9375ULL, 28178ULL, 450775ULL, 9780504ULL, 1795265022ULL }) {
            if (a % n == 0) continue;
            ll x = binpow(static_cast<ll>(a), d, n);
            if (x == 1 or x == n - 1) continue;
            bool composite = true;
            for (int r = 1; r < s; r++) {
                x = mul_mod(x, x, n);
                if (x == n - 1) {composite = false; break;}
            }
            if (composite) return false;
        }
        return true;
    }
    static ll pollard_rho(const ll n) {
        if (n % 2 == 0) return 2;
        if (n % 3 == 0) return 3;
        const ll c = uniform_int_distribution<ll>(1, n - 1)(rng);
        ll x = uniform_int_distribution<ll>(0, n - 1)(rng), y = x, d = 1;
        auto f = [&](const ll v) { return (mul_mod(v, v, n) + c) % n; };
        while (d == 1) {
            x = f(x), y = f(f(y)), d = gcd(llabs(x - y), n);
            if (d == n) return pollard_rho(n);
        }
        return d;
    }
    void factor(const ll n) {
        if (n == 1) return;
        if (is_prime(n)) factors[n]++;
        else { const ll d = pollard_rho(n); factor(d); factor(n / d); }
    }
};