//
// Created by heesang on 2/27/25.
//

#include <iostream>
#include <vector>
#include <chrono>

#include "Eigen/Sparse"
#include "Eigen/Dense"

#define MAX_N 10000



int main(int argc, char* argv[]) {
  Eigen::MatrixXd mat;

  // deterministic random sparse matrix generation; target sparsity = 0.1
  // https://stackoverflow.com/questions/30741884/sparse-random-matrix-with-eigen

  for (int N = 1; N <= MAX_N; N = N * 10) {


    std::cout << "N = " << N << std::endl;
    Eigen::MatrixXd mat1 = (Eigen::MatrixXd::Random(N, N).array() >= 0.8).cast<double>() * Eigen::MatrixXd::Random(N, N).array();
    Eigen::MatrixXd mat2 = (Eigen::MatrixXd::Random(N, N).array() >= 0.8).cast<double>() * Eigen::MatrixXd::Random(N, N).array();

    auto before = std::chrono::high_resolution_clock::now();
    Eigen::MatrixXd mat12 = mat1 * mat2;
    auto after = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double, std::milli> exec_time = after - before;
    std::cout << exec_time.count() << "ms" << std::endl;

    // dense

    // sparse column-major
    Eigen::SparseMatrix<double> spm1 = mat1.sparseView();
    Eigen::SparseMatrix<double> spm2 = mat2.sparseView();

    before = std::chrono::high_resolution_clock::now();
    Eigen::SparseMatrix<double> spm12 = spm1 * spm2;
    after = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double, std::milli> exec_time2 = after - before;
    std::cout << exec_time2.count() << "ms" << std::endl;

  }
  return 0;

}
