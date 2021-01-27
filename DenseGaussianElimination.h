#ifndef DENSEGAUSS_E_H
#define DENSEGAUSS_E_H

#include "GaussianElimination.h"
class DenseGaussE : public GaussE
{
private:
    /* data */
public:
    DenseGaussE(/* args */);
    ~DenseGaussE();
    std::vector<double> solve(Matrix& A, std::vector<double>& b);
};

#endif
