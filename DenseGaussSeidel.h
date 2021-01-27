#pragma once

#include "GaussSeidel.h"

class DenseGaussSeidel : public GaussSeidel
{
private:
    /* data */
public:
    DenseGaussSeidel(/* args */);
    ~DenseGaussSeidel();
    std::vector<double> solve(Matrix& A, std::vector<double>& b);
};