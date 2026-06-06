#include <bits/stdc++.h>
using namespace std;
using ll = long long int;

/*============================================================================================================
Binary Exponentiation, Permutation Powering & Fibonacci Utilities

Overview:
  This file collects several fast exponentiation–based techniques:
    • Modular exponentiation of integers
    • Repeated application of permutations using binary lifting
    • Fibonacci number computation using fast doubling and matrix exponentiation

  All algorithms run in logarithmic time with respect to the exponent
  and are suitable for competitive programming and number-theoretic tasks.

1. Modular Binary Exponentiation:
  - Computes (a^b) mod m efficiently using binary exponentiation
  - Reduces multiplication count from O(b) to O(log(b))
  - Handles large exponents safely under modulo arithmetic

2. Permutation Exponentiation:
  - Treats a permutation as a function and applies it k times
  - Uses the same binary exponentiation idea as modular power
  - apply_permutation:
      Applies a single permutation in O(n)
  - permute:
      Applies a permutation k times in O(n·log(k))

3. Fibonacci Numbers – Theory:
  - Definitions:
      F(0) = 0, F(1) = 1, F(n) = F(n - 1) + F(n - 2)
  - Key identities:
      • F(n−1)·F(n+1) − F(n)^2 = (−1)^n
      • F(n+k) = F(k)·F(n+1) + F(k−1)·F(n)
      • gcd(F(n), F(m)) = F(gcd(n, m))
      • Sum F(0..n) = F(n+2) − 1
      • Any natural number has a unique representation
          as a sum of non-consecutive Fibonacci numbers (Zeckendorf)

  - Generalized Fibonacci sum:
      For G(1) = x, G(2) = y:
        G(1..n) = x·F(n) + y·(F(n+1) − 1)

4. Fibonacci – Fast Doubling:
  - Computes {F(n), F(n+1)} recursively
  - Uses divide-and-conquer on binary representation of n
  - Highly efficient and numerically stable

5. Fibonacci – Matrix Exponentiation:
  - Uses the identity:
      | F(n+1) F(n)   | = base^n
      | F(n)   F(n−1) |

      where base = {{1,1},{1,0}}
  - Matrix exponentiation via binary powering

Conventions:
  • Permutations are 0-based
  • No modulo is applied in Fibonacci computations

Complexity Summary:
  • binpow:                O(log(b))
  • apply_permutation:     O(n)
  • permute:               O(n·log(k))
  • Fibonacci (doubling):  O(log(n))
  • Fibonacci (matrix):    O(log(n))
  • All algorithms use O(1) or O(n) extra space as appropriate

Applications:
  • Modular arithmetic and cryptography
  • Efficient simulation of repeated shuffles
  • Fibonacci-related number theory problems
  • Competitive programming and algorithm design

Notes:
  • binpow assumes mod > 0
  • Permutation exponentiation mutates the permutation vector
  • Fibonacci values grow exponentially and may overflow ll
============================================================================================================*/

ll binpow(ll a, ll b, const ll mod) {
    a %= mod;
    ll res = 1;
    while (b > 0) {if (b & 1) res = res * a % mod; a = a * a % mod, b >>= 1;}
    return res;
}

vector<int> apply_permutation(const vector<int> &sequence, const vector<int> &permutation) {
    vector<int> newSequence(sequence.size());
    for(int i = 0; i < sequence.size(); i++) newSequence[i] = sequence[permutation[i]];
    return newSequence;
}

vector<int> permute(vector<int> &sequence, vector<int> &permutation, ll k) {
    while (k > 0) {
        if (k & 1) sequence = apply_permutation(sequence, permutation);
        permutation = apply_permutation(permutation, permutation), k >>= 1;
    }
    return sequence;
}

pair<ll, ll> fib(const ll n) {
    if (n == 0) return {0, 1};
    auto [fst, snd] = fib(n >> 1);
    ll c = fst * (2 * snd - fst), d = fst * fst + snd * snd;
    if (n & 1) return {d, c + d};
    return {c, d};
}

struct matrix {
    ll mat[2][2];
    
    matrix friend operator *(const matrix &a, const matrix &b) {
        matrix c{};
        for (int i = 0; i < 2; i++) 
            for (int j = 0; j < 2; j++) {
                c.mat[i][j] = 0;
                for (int k = 0; k < 2; k++) c.mat[i][j] += a.mat[i][k] * b.mat[k][j];
            }
        return c;
    }
};

matrix matpow(matrix base, ll n) {
    matrix ans {{
        {1, 0},
        {0, 1}
    }};
    while (n > 0) { if(n & 1) ans = ans * base; base = base * base, n >>= 1; }
    return ans;
}

ll fib_mat(const ll n) {
    constexpr matrix base{{
        {1, 1},
        {1, 0}
    }};
    return matpow(base, n).mat[0][1];
}
