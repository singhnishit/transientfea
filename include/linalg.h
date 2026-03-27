#pragma once
#include <vector>
#include <string>

using Vec = std::vector<double>;
using Mat = std::vector<Vec>;

Mat mat_zeros(int n, int m);
Mat mat_eye(int n);
Vec vec_zeros(int n);
Vec vec_add(const Vec& a, const Vec& b);
Vec vec_sub(const Vec& a, const Vec& b);
Vec vec_scale(const Vec& v, double s);
Vec mat_mul_vec(const Mat& A, const Vec& v);
Mat mat_add(const Mat& A, const Mat& B);
Mat mat_scale(const Mat& A, double s);
Mat mat_add_scaled(const Mat& A, const Mat& B, double s); // A + s*B

// LU solve with partial pivoting: returns x such that A x = b
Vec lu_solve(Mat A, Vec b);

void print_vec(const Vec& v, const std::string& name = "");
void print_mat(const Mat& A, const std::string& name = "");
