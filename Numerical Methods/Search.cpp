#include <bits/stdc++.h>
using namespace std;
using ld = long double;
using ll = long long int;

/*============================================================================================================
Binary Search :

Description:
  • We use binary search to find a value in a sequence that is sorted
  • There are many ways to implement binary search
  • Another application of binary search is to find a place in a continuous function where f(l) <= f(r)
  • Also, most of the times, we do binary search on the answer of the problem

Notes:
  • In our implementation, we first define l and r out of the array bounds
  • We use the search interval [l, r)
  • At the end of the search:
    - l will be the index of the last element that is not greater than k (or -1 if there is no such element)
    - r will be the index of the first element larger than k (or n if there is no such element)
============================================================================================================*/

void binary_search(ll k, const vector<ll> &a) {
    int l = -1, r = a.size();
    while (r - l > 1) {
        int mid = l + (r - l) / 2;
        if (k < a[mid]) r = mid;
        else l = mid;
    }
}

/*============================================================================================================
Ternary Search :

Description:
  • Ternary Search (TS) is used to find the maximum or minimum of unimodal functions
  • A unimodal function is one that is first strictly increasing/decreasing and then strictly decreasing/increasing
  • The maximum or minimum point may lie over a continuous range of points

Types of Ternary Search:
  • Real numbers:
    - We define a fixed number of iterations (usually 200 to 300)
    - The search is performed on a continuous interval
  • Integer numbers:
    - Uses a different approach adapted for discrete values
    - This implementation finds the maximum value

Applications:
  • Finding optimal values in unimodal functions where derivative-based methods are difficult
  • Useful in optimization problems where the function shape is known or assumed unimodal

Notes:
  • The real-valued version repeatedly narrows the search interval based on function comparisons
  • The integer version narrows the search until the interval is small, then checks all candidates explicitly
============================================================================================================*/

ld f_real(ld x);

ll f_int(ll x);

ld ternary_search_real(ld l, ld r) {
    for (int i = 0; i < 200; i++) {
        ld m1 = l + (r - l) / 3, m2 = r - (r - l) / 3;
        if (f_real(m1) > f_real(m2)) r = m2;
        else l = m1;
    }
    return f_real(l);
}

ll ternary_search_int(ll l, ll r) {
    while (r - l > 4) {
        ll m1 = (l + r) / 2, m2 = (l + r) / 2 + 1;
        if (f_int(m1) > f_int(m2)) r = m2;
        else l = m1;
    }
    ll ans = -LLONG_MAX;
    for (ll i = l; i <= r; i++) ans = max(ans, f_int(i));
    return ans;
}
