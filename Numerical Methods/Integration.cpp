#include <bits/stdc++.h>
using namespace std;
typedef long double ld;

/*============================================================================================================
Numerical Integration: Composite Simpson's 1/3 Rule

Description:
  • Approximates the definite integral of a continuous function over the interval [a, b]
  • Utilizes Simpson's 1/3 Rule, a method that approximates the integrand by a quadratic polynomial
  • Suitable for smooth functions; accuracy improves with an increasing number of subintervals

Applications:
  • Computing areas under curves in physics, engineering, and economics
  • Estimating integrals in numerical simulations and data analysis

Notes:
  • The method divides the interval [a, b] into N subintervals (N must be even)
  • The formula is: 
    ∫[a,b] f(x) dx ≈ (h/3) * [f(a) + f(b) + 4 * Σ(f(x_i)) + 2 * Σ(f(x_j))]
    where x_i are odd-indexed points and x_j are even-indexed points within the interval
  • Convergence rate: O(h⁴), where h is the width of the subintervals

Order:
  • O(N), where N is the number of subintervals (note: N must be even)
============================================================================================================*/

ld func(ld x);

ld integration(ld a, ld b) {
    int N = 1000 * 1000;
    ld h = (b - a) / N;
    ld s = func(a) + func(b);
    for (int i = 1; i <= N - 1; i++) {
        ld x = a + h * i;
        s += func(x) * ((i & 1) ? 4 : 2);
    }
    s *= h / 3;
    return s;
}
