#include <bits/stdc++.h>
using namespace std;

typedef long double ld;
typedef long long int ll;
typedef pair<int, int> pii;
const int N = 1e6 + 10, ALPHA = 30, INF = 1e9 + 10, MOD = 1e9 + 7;

int f[N];
string s, t;
vector<int> pos;

void max_prefix_suffix() {
    int k = 0;
    f[0] = -1;
    f[1] = 0;
    for (int i = 2; i <= t.size(); i++) {
        while (k != -1 and t[k] != t[i - 1]) {
            k = f[k];
        }
        f[i] = ++k;
    }
}

void KMP() {
    max_prefix_suffix();
    int k = 0;
    for (int i = 0; i < s.size(); i++) {
        while (k != -1 and t[k] != s[i]) {
            k = f[k];
        }
        k++;
        if (k == t.size()) {
            k = f[k];
            pos.push_back(i + 2 - t.size());
        }
    }
}
