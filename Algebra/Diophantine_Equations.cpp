#include <bits/stdc++.h>
using namespace std;

/*--------------------------------------------------------------------------------------------------------
Find x and y such that: a.x + b.y = gcd(a, b)
Note that always: gcd(a, b) | a.x + b.y 
--------------------------------------------------------------------------------------------------------*/

int gcd(int a, int b, int &x, int &y) {
    if (b == 0) {
        x = 1, y = 0;
        return a;
    }
    int x1, y1;
    int g = gcd(b, a % b, x1, y1);
    x = y1;
    y = x1 - y1 * (a / b);
    return g;
}

/*--------------------------------------------------------------------------------------------------------
Solving: a.x + b.y = c
Note that: a != 0 and b != 0
--------------------------------------------------------------------------------------------------------*/

bool find_any_solution(int a, int b, int c, int &x, int &y, int &g) {
    g = gcd(abs(a), abs(b), x, y);
    if (c % g) return false;
    x *= c / g;
    y *= c / g;
    if (a < 0) x = -x;
    if (b < 0) y = -y;
    return true;    
}

/*--------------------------------------------------------------------------------------------------------
Solving "a.x + b.y = c" with range :)
return value : number of answers
minx <= x <= maxx
miny <= y <= maxy
After solve: lx <= x <= rx
x = lx + k.(b / g) and k >= 0
Order = Log(min(a, b))
For minimum x + y:
if a <= b => k = 0
if a > b => k = max
--------------------------------------------------------------------------------------------------------*/

void shift_solution(int &x, int &y, int a, int b, int cnt) {
    x += cnt * b;
    y -= cnt * a;
}

int find_all_solutions(int a, int b, int c, int minx, int maxx, int miny, int maxy) {
    int x, y, g;
    if (!find_any_solution(a, b, c, x, y, g)) return 0;
    a /= g, b /= g;
    int sign_a = a > 0 ? +1 : -1;
    int sign_b = b > 0 ? +1 : -1;

    shift_solution(x, y, a, b, (minx - x) / b);
    if (x < minx) shift_solution(x, y, a, b, sign_b);
    if (x > maxx) return 0;
    int lx1 = x;

    shift_solution(x, y, a, b, (maxx - x) / b);
    if (x > maxx) shift_solution(x, y, a, b, -sign_b);
    int rx1 = x;

    shift_solution(x, y, a, b, -(miny - y) / a);
    if (y < miny) shift_solution(x, y, a, b, -sign_a);
    if (y > maxy) return 0;
    int lx2 = x;

    shift_solution(x, y, a, b, -(maxy - y) / a);
    if (y > maxy) shift_solution(x, y, a, b, sign_a);
    int rx2 = x;

    if (lx2 > rx2) swap(lx2, rx2);
    int lx = max(lx1, lx2);
    int rx = min(rx1, rx2);

    if (lx > rx) return 0;
    return (rx - lx) / abs(b) + 1;
}
