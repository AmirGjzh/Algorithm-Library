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

int add(const int a, const int b, const int M) {
    const int res = (a % M + b % M) % M;
    return res < 0 ? res + M : res;
}

int sub(const int a, const int b, const int M) {
    const int res = (a % M - b % M) % M;
    return res < 0 ? res + M : res;
}

int mul(const int a, const int b, const int M) {
    const ll res = static_cast<ll>(a % M) * (b % M) % M;
    return res < 0 ? static_cast<int>(res) + M : static_cast<int>(res);
}

int ext(const int a, const int b, int &x, int &y) {
    if (b == 0) { x = 1; y = 0; return a; }
    int x1, y1;
    const int g = ext(b, a % b, x1, y1);
    x = y1;
    y = x1 - (a / b) * y1;
    return g;
}

int inverse(const int a, const int M) {
    int x, y;
    if (const int g = ext(a, M, x, y); g != 1) return -1;
    x = (x % M + M) % M;
    return x;
}

int div(const int a, const int b, const int M) {
    return static_cast<int>(static_cast<ll>(a % M + M) % M * inverse(b, M) % M);
}

int pow(int base, ll exp, const int M) {
    int res = 1;
    base = (base % M + M) % M;
    while (exp) {
        if (exp & 1) res = static_cast<int>(static_cast<ll>(res)) * base % M;
        base = static_cast<int>(static_cast<ll>(base)) * base % M;
        exp >>= 1;
    }
    return res;
}
