#include <bits/stdc++.h>
using namespace std;
using ll = long long int;

/*============================================================================================================
Bit Manipulation Utilities & Tricks

Overview:
  A compact collection of common bit manipulation operations and helper functions.
  Designed for competitive programming, low-level optimization, and bitmask-based algorithms.
  Covers single-bit operations, counting techniques, Gray code transformations,
  and submask enumeration patterns.

1. Basic Bit Operations (on integer n):
   1) Set the x-th bit:              n | (1 << x)
   2) Clear the x-th bit:            n & ~(1 << x)
   3) Flip the x-th bit:             n ^ (1 << x)
   4) Check the x-th bit:            (n >> x) & 1
   5) Check divisibility by 2^k:     (n & ((1 << k) - 1)) == 0
   6) Check if n is a power of two:  n > 0 && (n & (n - 1)) == 0
   7) Clear the rightmost set bit:   n & (n - 1)
   8) Clear all trailing ones:       n & (n + 1)
   9) Set the lowest cleared bit:    n | (n + 1)
   10) Extract the lowest set bit:   n & -n    (example: 00110100 → 00000100)

2. Counting Set Bits:
  - count_set_bits:
      Uses Brian Kernighan’s algorithm.
      Each iteration removes the lowest set bit.
  - count_set_bits_1_n:
      Counts total number of set bits in all numbers from 1 to n (inclusive)
      using a highest-set-bit decomposition approach.

3. Submask Enumeration:
  - finding_all_submasks:
      Iterates through all submasks of a given bitmask using:
        sub = (sub - 1) & mask
      Includes both mask itself and 0.

4. Gray Code:
  - gray:
      Converts a binary number to reflected Gray code.
      Adjacent numbers differ by exactly one bit.
  - gray_reverse:
      Converts a reflected Gray code back to its original binary value
      by cumulative XOR.

Conventions:
  • All bit positions are 0-based
  • Uses built-in intrinsics where appropriate

Complexity Summary:
  • count_set_bits:        O(number of set bits), worst-case O(log(n))
  • count_set_bits_1_n:    O(log(n))
  • Submask enumeration:   O(2^k) for k set bits in mask
  • Gray / Gray reverse:   O(log(n)) bit operations
  • All functions use O(1) extra space

Applications:
  • Bitmask DP and subset iteration
  • Low-level performance optimizations
  • Efficient counting problems
  • Gray codes in hardware design and algorithmic enumeration

Notes:
  • count_set_bits_1_n assumes n ≥ 0
  • Submask iteration pattern is standard and order-dependent
  • Gray codes minimize bit transitions, useful for state encoding
============================================================================================================*/

int count_set_bits(ll n) {
    int cnt = 0;
    while (n) n = n & n - 1, cnt++;
    return cnt;
}

ll count_set_bits_1_n(ll n) {
    if (n == 0) return 0;
    ll cnt = 0;
    while (n) {
        const ll x = 63 - __builtin_clzll(n);
        cnt += x << (x - 1);
        n -= 1LL << x;
        cnt += n + 1;
    }
    return cnt;
}

void finding_all_submasks(const ll mask) {
    for (ll sub = mask; ; sub = sub - 1 & mask) {
        // TODO
        if (sub == 0) break;
    }
}

ll gray(const ll n) {
    return n ^ n >> 1;
}

ll gray_reverse(ll g) {
    ll n = 0;
    for (; g; g >>= 1) n ^= g;
    return n;
}
