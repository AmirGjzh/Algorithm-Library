#include <bits/stdc++.h>
using namespace std;
using ll = long long int;

/*============================================================================================================
Euler’s Theorem & Number-Theory Utilities

Overview:
  This file contains classic arithmetic functions related to Euler’s theorem,
  divisor sums, totients, and general number-theoretic computations.

1. Euler’s Totient φ(n):
  - Definition: number of integers 1 ≤ i ≤ n with gcd(i, n) = 1
  - Formula: φ(n) = n · ∏_{p|n} (1 − 1/p) over all prime divisors p of n
  - Computation (trial division): O(√n)
  - Totient Sieve φ[1…n]: computes all φ(i) up to n
      • Initialize φ[i] = i
      • For each prime p ≤ n, update φ[j] -= φ[j]/p for multiples j of p
      • Complexity: O(n·log(log(n)))

2. Number of Divisors d(n):
  - If n = ∏ pᵉ, then d(n) = ∏ (e + 1)
  - Computation (trial division): O(√n)

3. Sum of Divisors σ(n):
  - If n = ∏ pᵉ, then σ(n) = ∏ (1 + p + p² + … + pᵉ)
  - Computation (trial division): O(√n)

4. Totient Summation over Divisors:
  - ∑_{d|n} φ(d) = n

5. Euler’s Theorem (for gcd(a, m) = 1):
  - a^φ(m) ≡ 1 (mod m)
  - Hence, for exponent n:
      a^n ≡ a^(n mod φ(m)) ≡ a^(φ(m) + (n mod φ(m))) (mod m)

6. Sum of gcd’s:
  - ∑_{i=1..n} gcd(i, n) = ∑_{d|n} d·φ(n/d)

7. Sum of lcm’s:
  - ∑_{i=1..n} lcm(i, n) = (n/2)·(1 + ∑_{d|n} d·φ(d))

Conventions:
  • Trial division used for factorization-based functions
  • Totient sieve assumes n ≥ 1

Complexity Summary:
  • phi(n)                      : O(√n)
  • phi_1_to_n(n)               : O(n·log(log(n)))
  • number_of_divisors(n)       : O(√n)
  • sum_of_divisors(n)          : O(√n)

Applications:
  • Modular arithmetic
  • Euler’s theorem problems
  • Divisor sums in number-theory algorithms
  • Competitive programming and problem-solving
============================================================================================================*/

ll phi(ll n) {
    ll result = n;
    for (ll i = 2; i * i <= n; i++) 
        if (n % i == 0) { while (n % i == 0) n /= i; result -= result / i; }
    if (n > 1) result -= result / n;
    return result;
}

vector<ll> phi_1_to_n(const ll n) {
    vector<ll> phi(n + 1);
    for (ll i = 0; i <= n; i++) phi[i] = i;
    for (ll i = 2; i <= n; i++) 
        if (phi[i] == i)  for (ll j = i; j <= n; j += i) phi[j] -= phi[j] / i;  
    return phi;                    
}

ll number_of_divisors(ll n) {
    ll total = 1;
    for (ll i = 2; i * i <= n; i++) 
        if (n % i == 0) {
            ll e = 0;
            do { e++, n /= i;
            } while (n % i == 0);
            total *= e + 1;
        }
    if (n > 1) total *= 2;
    return total;
}

ll sum_of_divisors(ll n) {
    ll total = 1;
    for (ll i = 2; i * i <= n; i++) 
        if (n % i == 0) {
            ll e = 0;
            do { e++, n /= i;
            } while (n % i == 0);
            ll sum = 0, pow = 1;
            do { sum += pow, pow *= i;
            } while (e-- > 0);
            total *= sum;
        }
    if (n > 1) total *= 1 + n;
    return total;
}
