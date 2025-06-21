#include <bits/stdc++.h>
using namespace std;
using u64 = uint64_t;
using u128 = __uint128_t;

/*============================================================================================================
Basic Primality Test (Trial Division)

Description:
  • Checks n ≤ √n by naive division — works in O(√n)
============================================================================================================*/

bool is_prime(int n) {
    for (int i = 2; i * i <= n; i++) 
        if (n % i == 0)
            return false;
    return n >= 2;
}

/*============================================================================================================
Probabilistic Primality Tests (Miller–Rabin / Fermat)

Description:
  • binpow: fast modular exponentiation
  • probably_prime_fermat: simple test based on Fermat's little theorem
  • miller_rabin1: standard probabilistic test using s, d decomposition
  • miller_rabin2: optimized deterministic version for 64-bit using fixed bases
============================================================================================================*/

struct ProbablyPrime {
    u64 binpow(u64 a, u64 b, u64 mod) {
        a %= mod;
        u64 res = 1;
        while (b > 0) {
            if (b & 1) res = (u128) res * a % mod;
            a = (u128) a * a % mod;
            b >>= 1;
        }
        return res;
    }

    bool probably_prime_fermat(int n, int iter = 5) {
        if (n < 4)
            return n == 2 or n == 3;
        for (int i = 0; i < iter; i++) {
            int a = 2 + rand() % (n - 3);
            if (binpow(a, n - 1, n) != 1)
                return false;
        }
        return true;
    }

    bool check_composite(u64 n, u64 a, u64 d, int s) {
        u64 x = binpow(a, d, n);
        if (x == 1 or x == n - 1)
            return false;
        for (int r = 1; r < s; r++) {
            x = (u128) x * x % n;
            if (x == n - 1)
                return false;
        }
        return true;
    }

    bool miller_rabin1(u64 n, int iter = 5) {
        if (n < 4)
            return n == 2 or n == 3;
        u64 s = 0, d = n - 1;
        while ((d & 1) == 0) 
            d >>= 1, s++;
        for (int i = 0; i < iter; i++) {
            int a = 2 + rand() % (n - 3);
            if (check_composite(n, a, d, s))
                return false;
        }
        return true;
    }

    bool miller_rabin2(u64 n) { 
        if (n < 2)
            return false;
        u64 r = 0, d = n - 1;
        while ((d & 1) == 0) 
            d >>= 1, r++;
        for (int a : {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37}) {
            if (n == a)
                return true;
            if (check_composite(n, a, d, r))
                return false;
        }
        return true;
    }
};

/*============================================================================================================
Sieve of Eratosthenes

Description:
  • Generates all primes up to n in O(n.log(log(n)))
  • Returns a boolean array of size n+1 where is_prime[i] is true if i is prime
============================================================================================================*/

vector<bool> sieve_of_eratosthenes(int n) {
    vector<bool> is_prime(n + 1, true);
    is_prime[0] = is_prime[1] = false;
    for (int i = 2; i * i <= n; i++) 
        if (is_prime[i]) 
            for (int j = i * i; j <= n; j += i)
                is_prime[j] = false;   
    return is_prime;                
}

/*============================================================================================================
Segmented Sieve

Description:
  • Finds primes in [L, R] by sieving using primes up to √R
  • Time: O((R – L).log(log(R)) + √R)
============================================================================================================*/

vector<char> segmented_sieve(int L, int R) {
    int lim = sqrt(R);
    vector<int> primes;
    vector<char> mark(lim + 1, false);
    for (int i = 2; i <= lim; i++) 
        if (!mark[i]) {
            primes.emplace_back(i);
            for (int j = i * i; j <= lim; j += i) 
                mark[j] = true;
        }
    vector<char> is_prime(R - L + 1, true);
    for (int i : primes) 
        for (int j = max(i * i, (L + i - 1) / i * i); j <= R; j += i)
            is_prime[j - L] = false;    
    if (L == 1) 
        is_prime[0] = false;
    return is_prime;            
}

/*============================================================================================================
Linear Sieve

Description:
  • Computes lp[i] = least prime divisor for all i ≤ n, in total O(n)
  • Also builds a list of primes
============================================================================================================*/

vector<int> linear_sieve(int n) {
    vector<int> primes;
    vector<int> lp(n + 1);
    for (int i = 2; i <= n; i++) {
        if (lp[i] == 0) {
            lp[i] = i;
            primes.push_back(i);
        }
        for (int j = 0; i * primes[j] <= n; j++) {
            lp[i * primes[j]] = primes[j];
            if (primes[j] == lp[i]) 
                break;
        }
    }
    return lp;
}

/*============================================================================================================
Trial Division Factorization

Description:
  • Extracts prime factors of n in O(√n)
  • Appends last prime > 1 if any remains
============================================================================================================*/

vector<int> trial_division(int n) {
    vector<int> factorization;
    for (int d = 2; d * d <= n; d++) 
        while (n % d == 0) {
            factorization.push_back(d);
            n /= d;
        }
    if (n > 1)
        factorization.push_back(n);
    return factorization;
}
