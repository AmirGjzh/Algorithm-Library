#include <bits/stdc++.h>
using namespace std;
using ll = long long int;

/*============================================================================================================
Description:
  Handy bit manipulation tricks and functions

  Basic bit operations:
    1 - Set the x-th bit:             n | (1 << x)
    2 - Clear the x-th bit:           n & ~(1 << x)
    3 - Flip the x-th bit:            n ^ (1 << x)
    4 - Check the x-th bit:           (n >> x) & 1
    5 - Check if n divisible by 2^k:  (n & (x - 1)) == 0
    6 - Check if n is a power of 2:   (n & (n - 1)) == 0
    7 - Clear the rightmost set bit:  n & (n - 1)
    8 - Clear all trailing ones:      n & (n + 1)
    9 - Set the last cleared bit:     n | (n + 1)
   10 - Extract the last set bit:     n & -n    (e.g., 00110100 → 00000100)

  Utility functions:
    - count_set_bits        : Brian Kernighan’s method (O(log(n)))
    - count_set_bits_1_n    : Count total set bits in [1..n] via bit-block approach
    - finding_all_submasks  : Iterate all submasks of a mask
    - gray                  : Convert n to reflected Gray code
    - gray_reverse          : Convert Gray code back to binary (inverse)
    - has_single_bit        : Checks if the number is a power of 2
    - bit_ceil / bit_floor  : Round up / down to the next power of two
    - rotl / rotr           : Rotate the bits in the number
    - countl_zero / countr_zero / countl_one / countr_one : Count the leading / trailing zeros / ones
    - popcount              : Count the number of set bits    

Time Complexity Summary:
  - count_set_bits:   O(log(n)) time, O(1) space
  - gray / gray_reverse: O(log(n)) bit operations

Applications:
  - Low-level programming, bitmask DP, competitive programming
  - Gray codes useful in hardware (error-minimizing encoders), algorithm design
============================================================================================================*/

int count_set_bits(ll n) {
    int cnt = 0;
    while (n) n = n & (n - 1), cnt++;
    return cnt;
}

ll count_set_bits_1_n(ll n) {
    if (n == 0) return 0;
    ll cnt = 0;
    while (n) {
        ll x = 63 - __builtin_clzll(n);
        cnt += x << (x - 1);
        n -= 1LL << x; 
        cnt += n + 1;
    }
    return cnt;
}

void finding_all_submasks(ll mask) {
    for (ll sub = mask; ; sub = (sub - 1) & mask) {
        // TODO
        if (sub == 0) break; // finished
    }
}

ll gray(ll n) {
    return n ^ (n >> 1);
}

ll gray_reverse(ll g) {
    ll n = 0;
    for (; g; g >>= 1) n ^= g;
    return n;    
}
