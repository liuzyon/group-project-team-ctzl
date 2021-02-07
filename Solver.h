#ifndef SOLVER_H
#define SOLVER_H

#include "Matrix.h"
#include "CSRMatrix.h"
#include <vector>

template <class T> class Solver
{
private:
    double get_error(Matrix<T>& A, T* b, T* x, int size, double error);
    T calMagnitude(T* vec, int size);
    void unitization(T* input, T* output, int size);
    T innerProduct(T* one, T* two, int size);

protected:
    /* data */
    // Matrix A;
    // std::vector<double> b;
public:
    Solver(/* args */);
    ~Solver();

//    virtual std::vector<double> solve(Matrix<T>& A, std::vector<double>& b) = 0;
    //**********iteration solver***************
    void dense_jacobi_solver(Matrix<T>& A, T* b, T* x, double tol = 1e-6);

    void dense_gauss_seidel_solver(Matrix<T>& A, T* b, T* x, double tol = 1e-6);

    void DenseGMRES(Matrix<T>& A, T* b, T* x);

    void dense_multigrid_solver(Matrix<T>& A, T* b, T* x);

    void sparse_jacobi_solver(CSRMatrix<T>& A, T* b, T* x, double tol = 1e-6);

    void sparse_gauss_seidel_solver(CSRMatrix<T>& A, T* b, T* x, double tol = 1e-6);

    // this one has not been finished yet
    void sparse_multigrid_solver(CSRMatrix<T>& A, T* b, T* x);
    //**********direct solver***************
    // Dense Gaussian Elimination solve
    void DenseGaussESolve(Matrix<T> A, T b[], T* x);

    void DenseGaussEPPSolve(Matrix<T> A, T b[], T* x);

    void DenseCholeskySolve(Matrix<T>& A, T* b, T* x);

    void DenseLUFactorisationSolve(Matrix<T>& A, T* b, T* x);

    // this one has not been finished yet
    void SparseLUFactorisationSolve(CSRMatrix<T>& A, T* b, T* x);
    // this one has not been finished yet
    void SparseCholeskySolve(CSRMatrix<T>& A, T* b, T* x);


    void printA();
    void printb();
};

#endif