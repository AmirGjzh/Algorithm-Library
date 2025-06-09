#include <bits/stdc++.h>
using namespace std;

/*--------------------------------------------------------------------------------------------------------
Binary Search :
We use binary search to find a value in a sequence that is sorted
There are many ways to do the implementation
Another application of BS is to find a place in a continuouse functions (f(l) <= f(r))
Also Most of the times, we do binary search on the answer of the problem

In our implementation we do these :
first we define l and r, out of the array bound
we use [l, r) 
at the end, l will be the index of the last element that is not greater than k (or -1 if there is no such element) and  
r will be the index of the first element larger than k  (or  n if there is no such element)
--------------------------------------------------------------------------------------------------------*/

void binary_search(int k, int l, int r, vector<int> &a) {
    l = -1, r = a.size();
    while (r - l > 1) {
        int mid = (l + r) >> 1;
        if (k < a[mid]) r = mid;
        else l = mid;
    }
}

/*--------------------------------------------------------------------------------------------------------
Ternary Search :
TS works just with unimodal functions
means they are first strictly inc/dec and then dec/inc
the max/min point can be some continues points

2 type of TS we have :
Real numbers => we define a constant and we do TS, that number of times (usually 200, 300)
Integer unmbers => we use a different method
We implement the maximum find
--------------------------------------------------------------------------------------------------------*/

double f_real(double x);

int f_int(int x);

double ternary_search_real(double l, double r) {
    for (int i = 0; i < 200; i++) {
        double m1 = l + (r - l) / 3;
        double m2 = r - (r - l) / 3;
        if (f_real(m1) > f_real(m2)) r = m2;
        else l = m1;
    }
    return f_real(l);
}

int ternary_search_int(int l, int r) {
    while (r - l > 4) {
        int m1 = (l + r) / 2;
        int m2 = (l + r) / 2 + 1;
        if (f_int(m1) > f_int(m2)) r = m2;
        else l = m1;
    }
    int ans = -1e9;
    for (int i = l; i <= r; i++) ans = max(ans, f_int(i));
    return ans;
}
