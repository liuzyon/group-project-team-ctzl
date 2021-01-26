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

    void printA();
    void printb();
};

#endif