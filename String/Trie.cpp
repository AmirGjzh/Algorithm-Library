#include <bits/stdc++.h>
using namespace std;

typedef long double ld;
typedef long long int ll;
typedef pair<int, int> pii;
const int N = 2e5 + 10, ALPHA = 30, INF = 1e9 + 10, MOD = 1e9 + 7;

int child[N][ALPHA], node_tag;
bool mark[N];

bool add_to_trie(string s) {
    bool new_node = false, flag = true;
    int cur = 0;
    for (int i = 0; i < s.size(); i++) {
        if (child[cur][s[i] - '0'] == 0) {
            child[cur][s[i] - '0'] = ++node_tag;
            new_node = true;
        }
        cur = child[cur][s[i] - '0'];
        if (mark[cur]) {
            flag = false;
        }
    }
    mark[cur] = true;
    if (!new_node) {
        flag = false;
    }
    return flag;
}
