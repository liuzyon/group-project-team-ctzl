#ifndef SOLVER_H
#define SOLVER_H

#include "Matrix.h"
#include <vector>

class Solver
{
protected:
    /* data */
    // Matrix A;
    // std::vector<double> b;
public:
    Solver(/* args */);
    ~Solver();

    virtual std::vector<double> solve(Matrix& A, std::vector<double>& b) = 0;

    // Dense Gaussian Elimination solve
    std::vector<double> DenseGaussESolve(Matrix& A, std::vector<double>& b);

    std::vector<double> DenseGaussEPPSolve(Matrix& A, std::vector<double>& b);

    std::vector<double> DenseGaussSeidelSolve(Matrix& A, std::vector<double>& b);

    std::vector<double> DenseJacobiSolve(Matrix& A, std::vector<double>& b);

    std::vector<double> DenseLUFactorisationSolve(Matrix& A, std::vector<double>& b);

    void printA();
    void printb();
};

#endif