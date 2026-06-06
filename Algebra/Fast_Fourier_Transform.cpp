#include <bits/stdc++.h>
using namespace std;
const double PI = acos(-1);
using ll = long long int;

/*============================================================================================================
Fast Fourier Transform (FFT) & Polynomial Multiplication

Overview:
  Implements polynomial multiplication using the Fast Fourier Transform.
  Polynomials are represented as coefficient vectors and multiplied via
  convolution in the frequency domain.

The implementation includes:
  • Recursive radix-2 FFT (educational reference)
  • Iterative radix-2 FFT with bit-reversal permutation (used in practice)
  • Polynomial multiplication using FFT and inverse FFT

1. FFT Strategy:
  - Input vectors are resized to the next power of two
  - FFT converts coefficient representation to point-value form
  - Pointwise multiplication is performed in O(n)
  - Inverse FFT reconstructs coefficients
  - Final values are normalized and rounded to integers

2. Implementations:
  - fft_recursive:
      Recursive Cooley–Tukey FFT implementation.
      Useful for understanding the algorithm structure.
  - fft:
      Iterative in-place FFT using bit-reversal permutation.
      More cache-friendly and efficient for large inputs.

3. Polynomial Multiplication:
  - multiply:
      Multiplies two integer polynomials using FFT-based convolution.
      Result size is rounded up to a power of two.
      Coefficients are rounded to the nearest integer after inverse FFT.

Conventions:
  • Polynomials are stored as vectors of coefficients in increasing degree order
  • Complex numbers use double precision
  • Angles are computed using 2·π / n roots of unity
  • No modular arithmetic, results rely on floating-point precision

Complexity Summary:
  • FFT / inverse FFT:     O(n·log(n))
  • Pointwise multiply:    O(n)
  • Polynomial multiply:   O(n·log(n))
  • Extra space:           O(n)

Applications:
  • Fast polynomial multiplication
  • Large integer multiplication via convolution
  • Signal processing
  • Computational algebra
  • Competitive programming

Notes:
  • Floating-point precision may introduce rounding errors for very large inputs
  • llround is used to convert final coefficients to integers
  • Input size should be kept within precision limits of double
============================================================================================================*/

void fft_recursive(vector<complex<double>> &a, const bool invert) {
    const int n = static_cast<int>(a.size());
    if (n == 1) return;
    vector<complex<double>> a0(n / 2), a1(n / 2);
    for (int i = 0; i * 2 < n; i++) a0[i] = a[2 * i], a1[i] = a[2 * i + 1];
    fft_recursive(a0, invert), fft_recursive(a1, invert);
    const double ang = 2 * PI / n * (invert ? -1 : 1);
    complex<double> w(1);
    const complex wn(cos(ang), sin(ang));
    for (int i = 0; i * 2 < n; i++) {
        a[i] = a0[i] + w * a1[i], a[i + n / 2] = a0[i] - w * a1[i];
        if (invert) a[i] /= 2, a[i + n / 2] /= 2; 
        w *= wn;
    }
}

void fft(vector<complex<double>> &a, const bool invert) {
    const int n = static_cast<int>(a.size());
    for (int i = 1, j = 0; i < n; i++) {
        int bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) swap(a[i], a[j]);    
    }
    for (int len = 2; len <= n; len <<= 1) {
        const double ang = 2 * PI / len * (invert ? -1 : 1);
        complex wlen(cos(ang), sin(ang));
        for (int i = 0; i < n; i += len) {
            complex<double> w(1);
            for (int j = 0; j < len / 2; j++) {
                complex<double> u = a[i + j], v = a[i + j + len / 2] * w;
                a[i + j] = u + v, a[i + j + len / 2] = u - v, w *= wlen;
            }
        }
    }
    if (invert) for (auto &x : a) x /= n;
}

vector<ll> multiply(const vector<ll> &a, const vector<ll> &b) {
    vector<complex<double>> fa(a.begin(), a.end()), fb(b.begin(), b.end());
    int n = 1;
    while (n < static_cast<int>(a.size() + b.size())) n <<= 1;
    fa.resize(n), fb.resize(n);
    fft(fa, false), fft(fb, false);
    for (int i = 0; i < n; i++) fa[i] *= fb[i];
    fft(fa, true);
    vector<ll> result(n);
    for (int i = 0; i < n; i++) result[i] = llround(fa[i].real());
    return result;
}
