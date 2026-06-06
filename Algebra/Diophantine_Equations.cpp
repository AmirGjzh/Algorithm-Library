#include <bits/stdc++.h>
using namespace std;
using ll = long long int;

/*============================================================================================================
Extended Euclidean Algorithm & Linear Diophantine Equations

Overview:
  This module provides a complete toolkit for working with linear Diophantine equations
  of the form:  a·x + b·y = c
  It includes the Extended Euclidean Algorithm, finding one valid solution,
  shifting solutions, and counting all integer solutions within given bounds.

1. Extended Euclidean Algorithm:
  - Computes gcd(a, b) and finds integers x and y such that:
      a·x + b·y = gcd(a, b)
  - Guarantees integer solutions for all integer inputs a, b
  - Forms the mathematical foundation for all subsequent functions

2. Solving a·x + b·y = c:
  - A solution exists iff c is divisible by gcd(a, b)
  - If a solution exists, one particular solution (x, y) is computed
  - Sign adjustments are applied to handle negative coefficients

3. General Solution Structure:
  - Given one solution (x₀, y₀) and g = gcd(a, b):
      x = x₀ + k·(b / g)
      y = y₀ − k·(a / g)
      for any integer k
  - Allows enumeration or shifting of solutions

4. Shifting Solutions:
  - shift_solution moves a known solution by k steps along the solution space
  - Used internally to align solutions with bounding constraints

5. Counting Bounded Solutions:
  - Counts all integer solutions (x, y) such that:
      minx ≤ x ≤ maxx
      miny ≤ y ≤ maxy
  - Efficiently computes the valid range of k values 
  - Returns the total number of valid (x, y) pairs

Conventions:
  • gcd(a, b) is always non-negative
  • Bounds are inclusive

Complexity Summary:
  • Extended GCD:            O(log(min(|a|, |b|)))
  • find_any_solution:       O(log(min(|a|, |b|)))
  • shift_solution:          O(1)
  • find_all_solutions:      O(log(min(|a|, |b|)))
  • All functions use O(1) extra space

Applications:
  • Solving linear Diophantine equations
  • Modular inverse computation
  • Cryptography (e.g., RSA, modular arithmetic)
  • Constraint satisfaction and counting integer solutions
  • Competitive programming and number theory problems

Notes:
  • If c % gcd(a, b) ≠ 0, no solution exists
  • Care must be taken with integer division and signs when shifting solutions
  • The solution count relies on precise boundary alignment
============================================================================================================*/

ll gcd(const ll a, const ll b, ll &x, ll &y) {
    if (b == 0) {x = 1, y = 0; return a;}
    ll x1, y1;
    const ll g = gcd(b, a % b, x1, y1);
    x = y1, y = x1 - y1 * (a / b);
    return g;
}

bool find_any_solution(const ll a, const ll b, const ll c, ll &x, ll &y, ll &g) {
    g = gcd(abs(a), abs(b), x, y);
    if (c % g) return false;
    x *= c / g, y *= c / g;
    if (a < 0) x = -x;
    if (b < 0) y = -y;
    return true;    
}

void shift_solution(ll &x, ll &y, const ll a, const ll b, const ll cnt) {
    x += cnt * b, y -= cnt * a;
}

ll find_all_solutions(ll a, ll b, const ll c, const ll minx, const ll maxx, const ll miny, const ll maxy) {
    ll x, y, g;
    if (!find_any_solution(a, b, c, x, y, g)) return 0;
    a /= g, b /= g;
    const ll sign_a = a > 0 ? +1 : -1;
    const ll sign_b = b > 0 ? +1 : -1;
    shift_solution(x, y, a, b, (minx - x) / b);
    if (x < minx) shift_solution(x, y, a, b, sign_b);
    if (x > maxx) return 0;
    const ll lx1 = x;
    shift_solution(x, y, a, b, (maxx - x) / b);
    if (x > maxx) shift_solution(x, y, a, b, -sign_b);
    const ll rx1 = x;
    shift_solution(x, y, a, b, -(miny - y) / a);
    if (y < miny) shift_solution(x, y, a, b, -sign_a);
    if (y > maxy) return 0;
    ll lx2 = x;
    shift_solution(x, y, a, b, -(maxy - y) / a);
    if (y > maxy) shift_solution(x, y, a, b, sign_a);
    ll rx2 = x;
    if (lx2 > rx2) swap(lx2, rx2);
    const ll lx = max(lx1, lx2);
    const ll rx = min(rx1, rx2);
    if (lx > rx) return 0;
    return (rx - lx) / abs(b) + 1;
}
