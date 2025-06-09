#include <bits/stdc++.h>
using namespace std;

/*--------------------------------------------------------------------------------------------------------
Tricks :
1 - set the x(th) bit : n | (1 << x)
2 - clear the x(th) bit : n & ~(1 << x)
3 - flip the x(th) bit : n ^ (1 << x)
4 - check the x(th) bit : (n >> x) & 1
5 - check if n is divisible by x = 2 ^ k : (n & (x - 1)) == 0
6 - check if n is a power of 2 : n & !(n & (n - 1))
7 - clear the rightmost set bit : n & (n - 1)
8 - clear all trailing ones : n & (n + 1)
9 - set the last cleared bit : n | (n + 1)
10 - extract the last set bit : n & -n   ex : 00110100 -> 00000100
Functions : 
has_single_bit : checks if the number is a power of 2
bit_ceil / bit_floor : round up / down to the next power of two
rotl / rotr : rotate the bits in the number
countl_zero / countr_zero / countl_one / countr_one : count the leading / trailing zeros / ones
popcount : count the number of set bits
--------------------------------------------------------------------------------------------------------*/

int count_set_bits(int n) {
    int cnt = 0;
    while (n) n = n & (n - 1), cnt++;
    return cnt;
}

int count_set_bits_1_n(int n) {
    int cnt = 0;
    while (n) {
        int x = __bit_width(n) - 1;
        cnt += x << (x - 1);
        n -= 1 << x;
        cnt += n + 1;
    }
    return cnt;
}

/*--------------------------------------------------------------------------------------------------------
Bitmask or mask : just a binary string :)
Submask of a mask : a mask that all its set bits are a subset of the original mask set bits
--------------------------------------------------------------------------------------------------------*/

void finding_all_submasks(int mask) {
    for (int sub = mask; ; sub = (sub - 1) & mask) {
        // TODO
        if (sub == 0) break; // finished
    }
}

/*--------------------------------------------------------------------------------------------------------
Convert n to its gray code and reverse
ex :
000 -> 000
001 -> 001
010 -> 011
011 -> 010
100 -> 110
101 -> 111
110 -> 101
111 -> 100
--------------------------------------------------------------------------------------------------------*/

int gray(int n) {
    return n ^ (n >> 1);
}

int gray_reverse(int g) {
    int n = 0;
    for (; g; g >>= 1)
        n ^= g;
    return n;    
}
