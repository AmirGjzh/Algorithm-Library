#include <bits/stdc++.h>
using namespace std;

/*--------------------------------------------------------------------------------------------------------
phi(i) = number of integers between 1 and  n  inclusive, which are coprime to n
ai | n ==> phi(a1) + phi(a2) + ... + phi(ak) = n 
Euler's theorem : a ^ phi(m) === 1  (mod m)     if m and a are coprime
a ^ n === a ^ (n % phi(m))  (mod m)     if m and a are coprime
a ^ n === a ^ (phi(m) + (n % phi(m)))  (mod m)
sum of (i : 1 => n) gcd(i, n) = sum of (d => d | n) (phi(n / d) * d)
sum of (i : 1 => n) lcm(i, n) = ((sum of (d => d | n) (phi(d) * d)) + 1) * n / 2
--------------------------------------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------------------------------------
Order = Sqrt(n)
Order = n.Log(Log(n))
--------------------------------------------------------------------------------------------------------*/

int phi(int n) {
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

void phi_1_to_n(int n) {
    vector<int> phi(n + 1);
    for (int i = 0; i <= n; i++)
        phi[i] = i;
    for (int i = 2; i <= n; i++) 
        if (phi[i] == i) 
            for (int j = i; j <= n; j += i)
                phi[j] -= phi[j] / i;          
}

/*--------------------------------------------------------------------------------------------------------
Order = Sqrt(n) (Both)
--------------------------------------------------------------------------------------------------------*/

int number_of_divisors(int n) {
    int total = 1;
    for (int i = 2; i * i <= n; i++) 
        if (n % i == 0) {
            int e = 0;
            do {
                e++;
                n /= i;
            } while (n % i == 0);
            total *= e + 1;
        }
    if (n > 1) total *= 2;
    return total;
}

int sum_of_divisors(int n) {
    int total = 1;
    for (int i = 2; i * i <= n; i++) 
        if (n % i == 0) {
            int e = 0;
            do {
                e++;
                n /= i;
            } while (n % i == 0);
            int sum = 0, pow = 1;
            do {
                sum += pow;
                pow *= i;
            } while (e-- > 0);
            total *= sum;
        }
    if (n > 1) total *= (1 + n);
    return total;
}
