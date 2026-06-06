#include <bits/stdc++.h>
using namespace std;
using ll = long long int;

/*============================================================================================================
Catalan Numbers

Description:
  • Catalan numbers C_n count many combinatorial objects.
    First few values (starting from C_0): 1, 1, 2, 5, 14, 42, 132, 429, 1430, ...
  • Two principal formulas / approaches:

  1) Recursive (convolution) DP:
     - C_0 = 1
     - For n >= 1: 
         C_n = sum_{k=0..n-1} C_k * C_{n-1-k}    (O(n²) DP)
     - Reliable for small/moderate n
     - Works for any MOD (no division issues)

  2) Analytical (closed form):
     - C_n = (1 / (n + 1)) * binom(2n, n)
     - Compute binomial using factorials, multiplicative formula, or big-integer libraries
     - Requires modular inverse if MOD is not prime
     - Fast for large n, but less flexible in CP if MOD not prime

Implementation Notes:
  • This function computes all Catalan numbers from C_0 to C_N modulo MOD
  • Suitable for competitive programming where many queries may be needed
  • Avoids division entirely, safe for arbitrary MOD
  • Time Complexity: O(n²)
  • Space Complexity: O(n)
  • For very large n (e.g., n > 1e6), analytical formula with modular arithmetic is preferred

Additional Notes:
  • The division by (n+1) in the exact formula is safe because binom(2n,n) is divisible by (n+1)
  • For modulo arithmetic: use modular inverse of (n+1) or precompute factorials/inverses

Applications:
  1) Number of correct bracket sequences with n pairs of parentheses
  2) Number of rooted full binary trees with n + 1 leaves (vertices not numbered)
     - A rooted binary tree is full if every vertex has either two children or no children
  3) Number of ways to completely parenthesize n + 1 factors
  4) Number of triangulations of a convex polygon with n + 2 sides
  5) Number of ways to connect 2n points on a circle with n disjoint chords (non-crossing)
  6) Number of non-isomorphic full binary trees with n internal nodes
     - Internal nodes = nodes with at least one child
  7) Number of monotonic lattice paths from (0,0) to (n,n) that never pass above y = x
     - Monotonic = steps (1,0) or (0,1)
  8) Number of permutations of length n that are stack-sortable (avoid pattern 3-1-2)
     - Equivalently: no indices i < j < k with a_k < a_i < a_j
  9) Number of non-crossing partitions of a set with n elements
 10) Number of ways to cover the ladder 1..n using n rectangles
     - i-th column has height i, known combinatorial interpretation
============================================================================================================*/

vector<ll> catalan(const int N, const int MOD) {
    vector<ll> C(N + 1, 0);
    C[0] = 1;
    for (int n = 1; n <= N; n++) {
        ll cur = 0;
        for (int k = 0; k <= n - 1; k++) {
            cur += C[k] * C[n - 1 - k];
            if (cur >= MOD) cur %= MOD; 
        }
        C[n] = cur % MOD;
    }
    return C;
}