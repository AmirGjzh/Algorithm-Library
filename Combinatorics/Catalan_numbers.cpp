#include <bits/stdc++.h>
using namespace std;
using ll = long long int;

/*============================================================================================================
Catalan Numbers

Description:
  • Catalan numbers C_n count many combinatorial objects
    First few values (starting from C_0): 1, 1, 2, 5, 14, 42, 132, 429, 1430, ...
  • Two principal formulas / approaches:
    1) Recursive (convolution): C_0 = C_1 = 1; for n >= 2:
      - C_n = sum_{k=0..n-1} C_k * C_{n-1-k}    (O(n^2) DP)
    2) Analytical (closed form): C_n = (1 / (n + 1)) * binom(2n, n)
      - compute binomial via factorials, multiplicative method, or big-integer libraries
  • Choose implementation by constraints:
    - Small n, many queries: precompute DP table O(N^2).
    - Large n but modulo prime and many queries: precompute factorials + inverse factorials (O(N) preprocess, O(1) query).
    - Exact C_n for moderate n where result fits in 64/128-bit: multiplicative binomial with gcd reduction then divide by (n+1).
    - Very large exact results: use big integers (boost::multiprecision::cpp_int).
  • Notes:
    - The division by (n+1) in exact integer formulas is safe because binom(2n,n) is divisible by (n+1)
    - For modulo arithmetic ensure you use modular inverse of (n+1) (or compute binom(2n,n) mod p and multiply by inv(n+1) mod p)
Applications:
  1) Number of correct bracket sequences consisting of n opening and n closing brackets           
  2) Number of rooted full binary trees with n + 1 leaves (vertices are not numbered)                   
    - A rooted binary tree is full if every vertex has either two children or no children
  3) Number of ways to completely parenthesize n + 1 factors                                    
  4) Number of triangulations of a convex polygon with n + 2 sides (partition into triangles)
  5) Number of ways to connect the 2n points on a circle to form n disjoint chords (non-crossing)
  6) Number of non-isomorphic full binary trees with n internal nodes (internal = node with at least one child)
  7) Number of monotonic lattice paths from (0,0) to (n,n) that do NOT pass above the main diagonal
    - Monotonic means steps (1,0) or (0,1); the path never goes above the diagonal y = x
  8) Number of permutations of length n that are stack-sortable (i.e., avoid pattern 3-1-2)           
    - Equivalently: no indices i<j<k with a_k < a_i < a_j
  9) Number of non-crossing partitions of a set of n elements
  10) Number of ways to cover the ladder 1..n using n rectangles (the i-th column has height i)    
    - Known combinatorial interpretation — see CP-Algorithms for the ladder-cover bijection
============================================================================================================*/

vector<ll> catalan(int N, ll MOD) {
    vector<ll> C(N + 1, 0);
    C[0] = 1;
    for (int n = 1; n <= N; n++) {
        __int128_t cur = 0;
        for (int k = 0; k <= n - 1; k++) {
            cur += (__int128_t) C[k] * C[n - 1 - k];
            if (cur >= MOD) cur %= MOD; 
        }
        C[n] = (ll) (cur % MOD);
    }
    return C;
}