#include <bits/stdc++.h>
using namespace std;
using ll = long long int;

/*============================================================================================================
Euler’s Theorem

1. Euler’s Totient φ(n)
   • Definition: count of integers 1 ≤ i ≤ n with gcd(i,n)=1
   • Formula: φ(n) = n·∏_{p|n} (1 − 1/p)
   • Computation (trial division): O(√n)

2. Totient Sieve φ[1…n]
   • Builds φ for all i≤n by initializing φ[i]=i and sieving primes:
       for p prime: for j=p,p*2,…: φ[j] -= φ[j]/p
   • Complexity: O(n.log(log(n)))

3. Number of Divisors d(n)
   • If n = ∏ pᵉ, then d(n) = ∏ (e+1)
   • Computation (trial division): O(√n)

4. Sum of Divisors σ(n)
   • If n = ∏ pᵉ, then σ(n) = ∏ (1 + p + p² + … + pᵉ)
   • Computation (trial division): O(√n)

5. Totient Summation over Divisors
   • ∑_{d|n} φ(d) = n

6. Euler’s Theorem (for gcd(a,m)=1)
   • a^φ(m) ≡ 1 (mod m)
   • Hence a^n ≡ a^(n mod φ(m)) ≡ a^(φ(m) + (n mod φ(m))) (mod m)

7. Sum of gcd’s:
   • ∑_{i=1..n} gcd(i,n) = ∑_{d|n} d·φ(n/d)

8. Sum of lcm’s:
   • ∑_{i=1..n} lcm(i,n) = (n/2)·(1 + ∑_{d|n} d·φ(d))
============================================================================================================*/

ll phi(ll n) {
    ll result = n;
    for (ll i = 2; i * i <= n; i++) 
        if (n % i == 0) {while (n % i == 0) n /= i; result -= result / i;}
    if (n > 1) result -= result / n;
    return result;
}

vector<ll> phi_1_to_n(ll n) {
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
    if (n > 1) total *= (1 + n);
    return total;
}
