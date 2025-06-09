#include <bits/stdc++.h>
using namespace std;

/*--------------------------------------------------------------------------------------------------------
Compute (a ^ b) % mod 
Order = Log(b)
--------------------------------------------------------------------------------------------------------*/

int binpow(int a, int b, int mod) {
    a %= mod;
    int res = 1;
    while (b > 0) {
        if (b & 1) res = (1LL * res * a) % mod;
        a = (1LL * a * a) % mod;
        b >>= 1;
    }
    return res;
}

/*-------------------------------------------------------------------------------------------------------- 
Applying a permutation to a sequence,  k times 
Order = Log(k) * sizeof(sequence)
--------------------------------------------------------------------------------------------------------*/

vector<int> apply_permutation(vector<int> sequence, vector<int> permutation) {
    vector<int> newSequence(sequence.size());
    for(int i = 0; i < sequence.size(); i++) 
        newSequence[i] = sequence[permutation[i]];
    return newSequence;
}

vector<int> permute(vector<int> sequence, vector<int> permutation, int k) {
    while (k > 0) {
        if (k & 1) 
            sequence = apply_permutation(sequence, permutation);
        permutation = apply_permutation(permutation, permutation);
        k >>= 1;
    }
    return sequence;
}

/*--------------------------------------------------------------------------------------------------------
Fibonacci 
Some rules:
F(0) = 0, F[1] = 1
F(n - 1).F(n + 1) - F(n) ^ 2 = (-1) ^ n
F(n + k) = F(k).F(n + 1) + F(k - 1).F(n)
F(n) | F(nk)
gcd(F(n), F(m)) = F(gcd(n, m))
F(0) + F(1) + ... + F(n) = F(n + 2) - 1

if G(1) = x and G(2) = y :
G(1) + G(2) + ... + G(n) = x * F(n) + y * (F(n + 1) - 1)

Any natural number N can be uniquely represented as a sum of Fibonacci numbers (solve with greedy)

{ F(n + 1), F(n) }  ==>  { 1, 1 } ^ n
{ F(n), F(n - 1) }  ==>  { 1, 0 }  

The sequence is periodic module M
--------------------------------------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------------------------------------
return {F(n), F(n + 1)}
Order = Log(n)
--------------------------------------------------------------------------------------------------------*/

pair<int, int> fib(int n) {
    if (n == 0) return {0, 1};
    auto p = fib(n >> 1);
    int c = p.first * (2 * p.second - p.first);
    int d = p.first * p.first + p.second * p.second;
    if (n & 1) return {d, c + d};
    return {c, d};
}

/*--------------------------------------------------------------------------------------------------------
return F(n)
Order = Log(n)
--------------------------------------------------------------------------------------------------------*/

struct matrix {
    int mat[2][2];
    matrix friend operator *(const matrix &a, const matrix &b) {
        matrix c;
        for (int i = 0; i < 2; i++) 
            for (int j = 0; j < 2; j++) {
                c.mat[i][j] = 0;
                for (int k = 0; k < 2; k++) 
                    c.mat[i][j] += a.mat[i][k] * b.mat[k][j];
            }
        return c;
    }
};

matrix matpow(matrix base, int n) {
    matrix ans {{
        {1, 0},
        {0, 1}
    }};
    while (n) {
        if(n & 1)
            ans = ans * base;
        base = base * base;
        n >>= 1;
    }
    return ans;
}

int fib_mat(int n) {
    matrix base{{
        {1, 1},
        {1, 0}
    }};
    return matpow(base, n).mat[0][1];
}
