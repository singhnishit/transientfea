#include "linalg.h"
#include <cmath>
#include <stdexcept>
#include <iostream>
#include <iomanip>

Mat mat_zeros(int n, int m) { return Mat(n, Vec(m, 0.0)); }
Mat mat_eye(int n) {
    Mat A = mat_zeros(n, n);
    for (int i = 0; i < n; i++) A[i][i] = 1.0;
    return A;
}
Vec vec_zeros(int n) { return Vec(n, 0.0); }

Vec vec_add(const Vec& a, const Vec& b) {
    Vec r(a.size());
    for (size_t i = 0; i < a.size(); i++) r[i] = a[i] + b[i];
    return r;
}
Vec vec_sub(const Vec& a, const Vec& b) {
    Vec r(a.size());
    for (size_t i = 0; i < a.size(); i++) r[i] = a[i] - b[i];
    return r;
}
Vec vec_scale(const Vec& v, double s) {
    Vec r(v.size());
    for (size_t i = 0; i < v.size(); i++) r[i] = v[i] * s;
    return r;
}
Vec mat_mul_vec(const Mat& A, const Vec& v) {
    int n = A.size();
    Vec r(n, 0.0);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < (int)v.size(); j++)
            r[i] += A[i][j] * v[j];
    return r;
}
Mat mat_add(const Mat& A, const Mat& B) {
    int n = A.size(), m = A[0].size();
    Mat C = mat_zeros(n, m);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            C[i][j] = A[i][j] + B[i][j];
    return C;
}
Mat mat_scale(const Mat& A, double s) {
    int n = A.size(), m = A[0].size();
    Mat C = mat_zeros(n, m);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            C[i][j] = A[i][j] * s;
    return C;
}
Mat mat_add_scaled(const Mat& A, const Mat& B, double s) {
    int n = A.size(), m = A[0].size();
    Mat C = mat_zeros(n, m);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            C[i][j] = A[i][j] + s * B[i][j];
    return C;
}

Vec lu_solve(Mat A, Vec b) {
    int n = A.size();
    std::vector<int> perm(n);
    for (int i = 0; i < n; i++) perm[i] = i;

    for (int k = 0; k < n; k++) {
        // partial pivot
        double maxval = std::abs(A[k][k]);
        int maxrow = k;
        for (int i = k+1; i < n; i++) {
            if (std::abs(A[i][k]) > maxval) { maxval = std::abs(A[i][k]); maxrow = i; }
        }
        std::swap(A[k], A[maxrow]);
        std::swap(b[k], b[maxrow]);
        std::swap(perm[k], perm[maxrow]);

        if (std::abs(A[k][k]) < 1e-14) continue;
        for (int i = k+1; i < n; i++) {
            double fac = A[i][k] / A[k][k];
            for (int j = k; j < n; j++) A[i][j] -= fac * A[k][j];
            b[i] -= fac * b[k];
        }
    }
    Vec x(n, 0.0);
    for (int i = n-1; i >= 0; i--) {
        x[i] = b[i];
        for (int j = i+1; j < n; j++) x[i] -= A[i][j] * x[j];
        if (std::abs(A[i][i]) > 1e-14) x[i] /= A[i][i];
    }
    return x;
}

void print_vec(const Vec& v, const std::string& name) {
    if (!name.empty()) std::cout << name << ": ";
    for (double x : v) std::cout << std::setw(12) << std::setprecision(4) << x << " ";
    std::cout << "\n";
}
void print_mat(const Mat& A, const std::string& name) {
    if (!name.empty()) std::cout << name << ":\n";
    for (auto& row : A) { for (double x : row) std::cout << std::setw(12) << std::setprecision(4) << x; std::cout << "\n"; }
}
