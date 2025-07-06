#include <bits/stdc++.h>
using namespace std;

/*============================================================================================================
Description:
  Compute (a ^ b) % mod efficiently using binary exponentiation

Time Complexity: O(log(b))

Applications:
  - Modular arithmetic in cryptography (e.g., RSA)
  - Competitive programming (modular inverses, large powers)
============================================================================================================*/

int binpow(int a, int b, int mod) {
    a %= mod;
    int res = 1;
    while (b > 0) {if (b & 1) res = (1LL * res * a) % mod; a = (1LL * a * a) % mod, b >>= 1;}
    return res;
}


/*============================================================================================================
Description:
  Apply a permutation multiple times efficiently using binary exponentiation on the permutation

Functions:
  - apply_permutation: Applies a single permutation
  - permute: Applies a permutation k times

Parameters:
  - sequence    : Vector of elements to permute
  - permutation : Permutation vector (0-based indices)
  - k           : Number of times to apply the permutation

Returns:
  New sequence after applying the permutation k times

Time Complexity:
  - apply_permutation : O(n)
  - permute           : O(n.log(k))

Applications:
  - Cyclic shifts and reordering patterns
  - Efficient simulation of repeated shuffles
============================================================================================================*/

vector<int> apply_permutation(vector<int> &sequence, vector<int> &permutation) {
    vector<int> newSequence(sequence.size());
    for(int i = 0; i < sequence.size(); i++) newSequence[i] = sequence[permutation[i]];
    return newSequence;
}

vector<int> permute(vector<int> &sequence, vector<int> &permutation, int k) {
    while (k > 0) {
        if (k & 1) sequence = apply_permutation(sequence, permutation);
        permutation = apply_permutation(permutation, permutation), k >>= 1;
    }
    return sequence;
}

/*============================================================================================================
Fibonacci Utilities

Rules & Identities:
  - F(0) = 0, F(1) = 1
  - F(n-1).F(n+1) - F(n)^2 = (-1)^n
  - F(n+k) = F(k).F(n+1) + F(k-1).F(n)
  - F(n) | F(n.k)
  - gcd(F(n), F(m)) = F(gcd(n, m))
  - Sum F(0 to n) = F(n + 2) - 1
  - Any natural number can be uniquely represented as a sum of non-consecutive Fibonacci numbers

Additional Series Identity:
  If a sequence G is defined such that:
    - G(1) = x
    - G(2) = y

  Then the sum G(1) + G(2) + ... + G(n) is given by:
    G(1 to n) = x.F(n) + y.(F(n + 1) - 1)  
============================================================================================================*/

/*============================================================================================================
Description:
  Compute Fibonacci numbers using the fast doubling method

Returns:
  Pair {F(n), F(n+1)}

Time Complexity: O(log(n))
============================================================================================================*/

pair<int, int> fib(int n) {
    if (n == 0) return {0, 1};
    auto p = fib(n >> 1);
    int c = p.first * (2 * p.second - p.first), d = p.first * p.first + p.second * p.second;
    if (n & 1) return {d, c + d};
    return {c, d};
}

/*===========================================================================================
Description:
  Compute the n-th Fibonacci number using matrix exponentiation

Matrix Form:
  | F(n+1) F(n)   |
  | F(n)   F(n-1) | = (base)^n, where base = {{1, 1}, {1, 0}}

Time Complexity: O(log(n))
===========================================================================================*/

struct matrix {
    int mat[2][2];
    matrix friend operator *(const matrix &a, const matrix &b) {
        matrix c;
        for (int i = 0; i < 2; i++) 
            for (int j = 0; j < 2; j++) {
                c.mat[i][j] = 0;
                for (int k = 0; k < 2; k++) c.mat[i][j] += a.mat[i][k] * b.mat[k][j];
            }
        return c;
    }
};

matrix matpow(matrix base, int n) {
    matrix ans {{
        {1, 0},
        {0, 1}
    }};
    while (n) {if(n & 1) ans = ans * base; base = base * base, n >>= 1;}
    return ans;
}

int fib_mat(int n) {
    matrix base{{
        {1, 1},
        {1, 0}
    }};
    return matpow(base, n).mat[0][1];
}
