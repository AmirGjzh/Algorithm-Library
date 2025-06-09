#include <bits/stdc++.h>
using namespace std;
const double PI = acos(-1);

/*--------------------------------------------------------------------------------------------------------
Multiplying to polynomials :
we store all coefficients of all powers in vector for each polynomial
the number of coefficients are power of 2
Order = n.Log(n)
--------------------------------------------------------------------------------------------------------*/

void fft_recursive(vector<complex<double>> &a, bool invert) {
    int n = a.size();
    if (n == 1) return;
    vector<complex<double>> a0(n / 2), a1(n / 2);
    for (int i = 0; i * 2 < n; i++) {
        a0[i] = a[2 * i];
        a1[i] = a[2 * i + 1];
    }
    fft_recursive(a0, invert);
    fft_recursive(a1, invert);
    double ang = 2 * PI / n * (invert ? -1 : 1);
    complex<double> w(1), wn(cos(ang), sin(ang));
    for (int i = 0; i * 2 < n; i++) {
        a[i] = a0[i] + w * a1[i];
        a[i + n / 2] = a0[i] - w * a1[i];
        if (invert) {
            a[i] /= 2;
            a[i + n / 2] /= 2;
        } 
        w *= wn;
    }
}

void fft(vector<complex<double>> &a, bool invert) {
    int n = a.size();
    for (int i = 1, j = 0; i < n; i++) {
        int bit = n >> 1;
        for (; j & bit; bit >>= 1)
            j ^= bit;
        j ^= bit;
        if (i < j)
            swap(a[i], a[j]);    
    }
    for (int len = 2; len <= n; len <<= 1) {
        double ang = 2 * PI / len * (invert ? -1 : 1);
        complex<double> wlen(cos(ang), sin(ang));
        for (int i = 0; i < n; i += len) {
            complex<double> w(1);
            for (int j = 0; j < len / 2; j++) {
                complex<double> u = a[i + j], v = a[i + j + len / 2] * w;
                a[i + j] = u + v;
                a[i + j + len / 2] = u - v;
                w *= wlen;
            }
        }
    }
    if (invert)
        for (auto &x : a)
            x /= n;
}

vector<int> multiply(vector<int> const &a, vector<int> const &b) {
    vector<complex<double>> fa(a.begin(), a.end()), fb(b.begin(), b.end());
    int n = 1;
    while (n < a.size() + b.size()) 
        n <<= 1;
    fa.resize(n);
    fb.resize(n);
    fft(fa, false);
    fft(fb, false);
    for (int i = 0; i < n; i++)
        fa[i] *= fb[i];
    fft(fa, true);
    vector<int> result(n);
    for (int i = 0; i < n; i++)
        result[i] = round(fa[i].real());
    return result;
}
