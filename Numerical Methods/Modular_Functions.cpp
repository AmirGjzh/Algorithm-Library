#include <bits/stdc++.h>
using namespace std;
using ll = long long int;

/*============================================================================================================
Modular Arithmetic:

Description:
  • Provides standalone functions for modular addition, subtraction, multiplication, division (via modular inverse), and exponentiation
  • Ensures results are always in the range [0, M-1], handling negative inputs gracefully
  • Uses the extended Euclidean algorithm to compute modular inverses for general (possibly composite) M
  • Includes fast binary exponentiation for powering under a modulus

Notes:
  • add(a, b, M): returns (a + b) mod M with non-negative output
  • sub(a, b, M): returns (a - b) mod M with non-negative output
  • mul(a, b, M): returns (a * b) mod M with non-negative output
  • div(a, b, M): returns a * inv(b) mod M, inv(b) computed via ext to handle gcd = 1
  • pow(base, exp, M): binary exponentiation to compute base^exp mod M efficiently
============================================================================================================*/

int add(int a, int b, int M) {
    int res = (a % M + b % M) % M;
    return res < 0 ? res + M : res;
}

int sub(int a, int b, int M) {
    int res = (a % M - b % M) % M;
    return res < 0 ? res + M : res;
}

int mul(int a, int b, int M) {
    ll res = (ll) (a % M) * (b % M) % M;
    return res < 0 ? res + M : res;
}

int div(int a, int b, int M) {
    return (int) ((ll) (a % M + M) % M * inverse(b, M) % M);
}

int ext(int a, int b, int &x, int &y) {
    if (b == 0) { x = 1; y = 0; return a; }
    int x1, y1;
    int g = ext(b, a % b, x1, y1);
    x = y1;
    y = x1 - (a / b) * y1;
    return g;
}

int inverse(int a, int M) {
    int x, y;
    int g = ext(a, M, x, y);
    if (g != 1) return -1; 
    x = (x % M + M) % M;
    return x;
}

int pow(int base, ll exp, int M) {
    int res = 1;
    base = (base % M + M) % M;
    while (exp) {
        if (exp & 1) res = (ll) res * base % M;
        base = (ll) base * base % M;
        exp >>= 1;
    }
    return int(res);
}
