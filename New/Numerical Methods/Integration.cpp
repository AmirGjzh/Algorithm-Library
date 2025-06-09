#include <bits/stdc++.h>
using namespace std;
typedef long double ld;

/*--------------------------------------------------------------------------------------------------------
Calculating integration of a function, in interval [a, b]
Just it :)
--------------------------------------------------------------------------------------------------------*/

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
