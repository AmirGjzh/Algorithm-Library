#include <bits/stdc++.h>
using namespace std;

/*============================================================================================================
Modular Inverse (single a, prime m)

Find x such that a·x ≡ 1 (mod m), Exists iff gcd(a,m)=1

Methods:
  • Euler’s: x ≡ a^(φ(m)−1)  (mod m)  when m is prime or gcd(a,m)=1
  • Extended GCD: solve a·x + m·y = 1

This recursive implementation works in O(log(m)) time for prime m:
  inv(a) = m − ⌊m/a⌋·inv(m mod a)  (mod m)  
============================================================================================================*/

int modular_inverse(int a, int m) {
    a %= m;
    return a <= 1 ? a : m - (m / a) * modular_inverse(m % a, m) % m;
}

/*============================================================================================================
Modular Inverses for 1 … m−1 (prime m)

Compute inv[i] = i^(m−2) mod m for all 1 ≤ i < m in O(m) total:
  inv[1] = 1
  inv[i] = m − ⌊m/i⌋·inv[m mod i]  (mod m)

Order: O(m)
============================================================================================================*/

vector<int> modular_inverse_1_m(int m) {
    vector<int> inv(m);
    inv[1] = 1;
    for (int i = 2; i < m; i++) 
        inv[i] = m - (m / i) * inv[m % i] % m;
    return inv;    
}

/*============================================================================================================
Chinese Remainder Theorem (pairwise coprime moduli)

Solve system of congruences:
  x ≡ a_i (mod m_i),    for i = 1 … k

Conditions:
  • Moduli m_i must be pairwise coprime
    – If not, you must first merge/consolidate congruences:
      • For any i≠j, let g = gcd(m_i, m_j). You need g | (a_i – a_j) or there is no solution
      • When g > 1, you can combine those two congruences into one with modulus lcm(m_i, m_j)
  • Input format: vector of pairs `{a_i, m_i}`

Algorithm:
  1. Compute M = ∏ m_i.
  2. For each (a_i, m_i):
      M_i = M / m_i
      N_i = modular_inverse(M_i mod m_i, m_i)
      term = a_i * M_i * N_i
      Accumulate `solution = (solution + term) mod M`
  3. Return `solution`

Complexity:  
  • O(k + log(M))
============================================================================================================*/

int chinese_remainder_theorem(vector<pair<int, int>> congruences) {
    int M = 1, solution = 0;
    for (auto &c : congruences) 
        M *= c.second;
    for (auto &c : congruences) {
        int a_i = c.first;
        int M_i = M / c.second;
        int N_i = modular_inverse(M_i, c.second);
        solution = (solution + a_i * M_i % M * N_i) % M;
    }
    return solution;
}

/*============================================================================================================
Discrete Logarithm (Baby‑Step Giant‑Step)

Find smallest x ≥ 0 such that a^x ≡ b (mod m), if it exists
Order: O(√m.log(m)) due to maps and modular multiplications
============================================================================================================*/

int baby_step_giant_step(int a, int b, int m) {
    a %= m, b %= m;
    if (a == 0) 
        return b == 0 ? 1 : -1;
    int k = 1, add = 0, g;
    while ((g = __gcd(a, m)) > 1) {
        if (b == k)
            return add;
        if (b % g)
            return -1;
        b /= g, m /= g, add++;
        k = (k * a / g) % m;        
    }
    int n = sqrt(m) + 1;
    int an = 1;
    for (int i = 0; i < n; i++)
        an = an * a % m;
    unordered_map<int, int> vals;
    for (int q = 0, cur = b; q <= n; q++) {
        vals[cur] = q;
        cur = cur * a % m;
    }
    for (int p = 1, cur = k; p <= n; p++) {
        cur = cur * an % m;
        if (vals.count(cur)) {
            int ans = n * p - vals[cur] + add;
            return ans;
        }
    }
    return -1;
}

/*============================================================================================================
Primitive Root modulo p (prime)

Description:
  • A number g is a primitive root modulo p if for every a with gcd(a,p)=1 there exists k such that:  
    g^k ≡ a  (mod p)
  • Equivalently, the powers of g (mod p) cycle through all nonzero residues before repeating

Existence:
  • A primitive root modulo n exists if and only if n is one of:
      1. 1, 2, or 4
      2. p^k for an odd prime p
      3. 2·p^k for an odd prime p
  • In particular, for a prime p there always exists a primitive root

Count:
  • The number of primitive roots modulo n (when they exist) is φ(φ(n))

Key Fact:
  • g is a primitive root mod n ⇔ the smallest k>0 with g^k ≡ 1 (mod n) is exactly φ(n)
  • This function finds the smallest primitive root for a prime p

Order: O((number of distinct factors of φ) · log p), typically fast
============================================================================================================*/

int binpow(int a, int b, int mod) {
    a %= mod;
    int res = 1;
    while (b > 0) {
        if (b & 1) res = res * a % mod;
        a = a * a % mod;
        b >>= 1;
    }
    return res;
}

int Phi(int n) {
    int result = n;
    for (int i = 2; i * i <= n; i++) {
        if (n % i == 0) {
            while (n % i == 0)
                n /= i;
            result -= result / i;
        }
    }
    if (n > 1)
        result -= result / n;
    return result;
}

int find_primitive_root(int p) {
    vector<int> fact;
    int phi = Phi(p), n = phi;
    for (int i = 2; i * i <= n; i++)
        if (n % i == 0) {
            fact.push_back(i);
            while (n % i == 0)
                n /= i;
        }
    if (n > 1)
        fact.push_back(n);
    for (int res = 2; res <= p; res++) {
        if (__gcd(res, p) != 1) continue;
        bool ok = true;
        for (size_t i = 0; i < fact.size() && ok; i++)
            ok &= (binpow(res, phi / fact[i], p) != 1);
        if (ok) return res;
    }
    return -1;
}

/*============================================================================================================
Discrete k‑th Root modulo prime n

Solve x^k ≡ a (mod n) for all x ∈ [0…n−1], n must be prime
Complexity: around O(√n·log n)
============================================================================================================*/

vector<int> find_all_discrete_roots(int n, int k, int a) {
    if (a == 0) return {0};
    int g = find_primitive_root(n);
    int sq = sqrt(n) + 1;
    vector<pair<int, int>> dec(sq);
    for (int i = 1; i <= sq; i++) 
        dec[i - 1] = {binpow(g, i * sq * k % (n - 1), n), i};
    sort(dec.begin(), dec.end());
    int any_ans = -1;
    for (int i = 0; i < sq; i++) {
        int my = binpow(g, i * k % (n - 1), n) * a % n;
        auto it = lower_bound(dec.begin(), dec.end(), make_pair(my, 0));
        if (it != dec.end() and it->first == my) {
            any_ans = it->second * sq - i;
            break;
        }
    }    
    if (any_ans == -1) return {};
    int delta = (n - 1) / __gcd(k, n - 1);
    vector<int> ans;
    for (int cur = any_ans % delta; cur < n - 1; cur += delta)
        ans.push_back(binpow(g, cur, n));
    sort(ans.begin(), ans.end());
    return ans;    
}
