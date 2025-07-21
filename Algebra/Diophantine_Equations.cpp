#include <bits/stdc++.h>
using namespace std;
using ll = long long int;

/*============================================================================================================
Description:
  Extended Euclidean Algorithm
  Finds integers x and y such that: a.x + b.y = gcd(a, b)

Mathematical Identity:
  - The equation a.x + b.y = gcd(a, b) always has integer solutions
  - Useful in solving Diophantine equations and computing modular inverses

Time Complexity: O(log(min(a, b)))

Applications:
  - Solving linear Diophantine equations
  - Finding modular inverse of a modulo b
  - Cryptography algorithms (e.g., RSA)
============================================================================================================*/

ll gcd(ll a, ll b, ll &x, ll &y) {
    if (b == 0) {x = 1, y = 0; return a;}
    ll x1, y1, g = gcd(b, a % b, x1, y1);
    x = y1, y = x1 - y1 * (a / b);
    return g;
}

/*============================================================================================================
Description:
  Solves the Diophantine equation a.x + b.y = c for any integers a, b, c (a ≠ 0, b ≠ 0)

Returns:
  - true if a solution exists, false otherwise
  - Outputs one valid solution (x, y) and gcd(a, b)

Time Complexity: O(log(min(a, b)))

Applications:
  - Finding any solution to linear Diophantine equations
  - Useful when only one solution is required, or to generate general solutions
============================================================================================================*/

bool find_any_solution(ll a, ll b, ll c, ll &x, ll &y, ll &g) {
    g = gcd(abs(a), abs(b), x, y);
    if (c % g) return false;
    x *= c / g, y *= c / g;
    if (a < 0) x = -x;
    if (b < 0) y = -y;
    return true;    
}

/*============================================================================================================
Description:
  Shifts a known solution (x, y) of a.x + b.y = c to another solution by k steps
  General solution form: x = lx + k.(b/g), y = ly - k.(a/g)

Time Complexity: O(1)

Applications:
  - Used in adjusting the solution to fall within a specific range
============================================================================================================*/


/*============================================================================================================
Description:
  Finds all integer solutions to the equation a.x + b.y = c that satisfy:
    minx <= x <= maxx
    miny <= y <= maxy

Returns:
  - The number of valid (x, y) integer pairs within the bounds

Time Complexity: O(log(min(a, b)))

Applications:
  - Counting constrained integer solutions to Diophantine equations
  - Resource allocation, cryptographic constraint solving
============================================================================================================*/

void shift_solution(ll &x, ll &y, ll a, ll b, ll cnt) {
    x += cnt * b, y -= cnt * a;
}

ll find_all_solutions(ll a, ll b, ll c, ll minx, ll maxx, ll miny, ll maxy) {
    ll x, y, g;
    if (!find_any_solution(a, b, c, x, y, g)) return 0;
    a /= g, b /= g;
    ll sign_a = a > 0 ? +1 : -1, sign_b = b > 0 ? +1 : -1;
    shift_solution(x, y, a, b, (minx - x) / b);
    if (x < minx) shift_solution(x, y, a, b, sign_b);
    if (x > maxx) return 0;
    ll lx1 = x;
    shift_solution(x, y, a, b, (maxx - x) / b);
    if (x > maxx) shift_solution(x, y, a, b, -sign_b);
    ll rx1 = x;
    shift_solution(x, y, a, b, -(miny - y) / a);
    if (y < miny) shift_solution(x, y, a, b, -sign_a);
    if (y > maxy) return 0;
    ll lx2 = x;
    shift_solution(x, y, a, b, -(maxy - y) / a);
    if (y > maxy) shift_solution(x, y, a, b, sign_a);
    ll rx2 = x;
    if (lx2 > rx2) swap(lx2, rx2);
    ll lx = max(lx1, lx2), rx = min(rx1, rx2);
    if (lx > rx) return 0;
    return (rx - lx) / abs(b) + 1;
}
