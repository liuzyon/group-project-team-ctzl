#ifndef SOLVER_H
#define SOLVER_H

#include "Matrix.h"
#include <vector>

template <class T> class Solver
{
protected:
    /* data */
    // Matrix A;
    // std::vector<double> b;
public:
    Solver(/* args */);
    ~Solver();

    virtual std::vector<double> solve(Matrix<T>& A, std::vector<double>& b) = 0;

    // Dense Gaussian Elimination solve
    std::vector<double> DenseGaussESolve(Matrix<T>& A, std::vector<double>& b);

    std::vector<double> DenseGaussEPPSolve(Matrix<T>& A, std::vector<double>& b);

    std::vector<double> DenseGaussSeidelSolve(Matrix<T>& A, std::vector<double>& b);

    std::vector<double> DenseJacobiSolve(Matrix<T>& A, std::vector<double>& b);

    std::vector<double> DenseLUFactorisationSolve(Matrix<T>& A, std::vector<double>& b);

    void printA();
    void printb();
};

#endif