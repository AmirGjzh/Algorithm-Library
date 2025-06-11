#include <bits/stdc++.h>
using namespace std;

const int N = 2e5 + 10, INF = 1e9 + 10, MOD = 1e9 + 10;

int _pow(int x, int e) {
    if (e == 0) {
        return 1;
    }
    return _pow(x * x % MOD, e / 2) * (e & 2 ? x  : 1);
}

int rev(int x) {
    return _pow(x, MOD - 2);
}