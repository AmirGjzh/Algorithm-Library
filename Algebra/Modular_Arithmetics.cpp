#include <bits/stdc++.h>
using namespace std;

/*--------------------------------------------------------------------------------------------------------
Modular Inverse 'a' :
such x that : a.x === 1 (mod m)
it exists if and only if gcd(a, m) == 1
FIRST WAY =>  x === a ^ (phi(m) - 1) (mod m)   when a and m are coprime
SECOND WAY =>  solve this : a.x + m.y = 1
A function that calculate this for PRIME m (very fast) :
--------------------------------------------------------------------------------------------------------*/

int modular_inverse(int a, int m) {
    a %= m;
    return a <= 1 ? a : m - (m / a) * modular_inverse(m % a, m) % m;
}

/*--------------------------------------------------------------------------------------------------------
Finding modular inverse of all numbers in [1, ..., m - 1] (m is prime)
Order = m
--------------------------------------------------------------------------------------------------------*/

void modular_inverse_1_m(int m) {
    vector<int> inv(m);
    inv[1] = 1;
    for (int i = 2; i < m; i++) 
        inv[i] = m - (m / i) * inv[m % i] % m;
}

/*--------------------------------------------------------------------------------------------------------
For solving these type of equations :
a.x === b (mod m)
we can solve this : a.x + m.y = b
--------------------------------------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------------------------------------
Chinese Remainder Theorem :
a === a1 (mod m1)
a === a2 (mod m2)
    .
    .
    .
a === ak (mod mk)
mi are coprime (if not, we must make it and check the conflicts)
what is a (mod (m1)(m2)...(mk))?   
pair<int, int> => {ai, mi} 
--------------------------------------------------------------------------------------------------------*/

int chinese_remainder_theorem(vector<pair<int, int>> congruences) {
    int M = 1, solution = 0;
    for (auto c : congruences) 
        M *= c.second;
    for (auto c : congruences) {
        int a_i = c.first;
        int M_i = M / c.second;
        int N_i = modular_inverse(M_i, c.second);
        solution = (solution + a_i * M_i % M * N_i) % M;
    }
    return solution;
}

/*--------------------------------------------------------------------------------------------------------
Discrete Logarithm :
such x that : a ^ x === b (mod m)     (does not always exist)
we can find the minimum x with this function :
Order = Sqrt(m)
--------------------------------------------------------------------------------------------------------*/

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

/*--------------------------------------------------------------------------------------------------------
Primitive Root :
A number g is a primitive root module n, if every number coprime to n is congruent to a power of  g modulo  n
for every a that gcd(a, n) = 1, there is a k that : g ^ k === a (mod n)
Primitive root modulo  n  exists if and only if :
1 - n is 1, 2, 4 or
2 - n is p ^ k (p is prime and odd)
3 - n is 2 * (p ^ k) (p is prime and odd)
the number of primitive roots module n, is phi(phi(n))
ONE Fact :
if g is a primitive root module n ==> the smallest k that g ^ k === 1 (mod n) is equal to phi(n)
and the reverse is true too
this function finds a primitive root module p (the smallest)
Order = ans.Log(phi(n)).Log(n)
--------------------------------------------------------------------------------------------------------*/

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

/*--------------------------------------------------------------------------------------------------------
Discrete Root :
Finding all x such that : x ^ k === a (mod n) (n is prime)
note that 0 <= x <= n - 1
it may have no answer
Order = don't know :)
--------------------------------------------------------------------------------------------------------*/

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
