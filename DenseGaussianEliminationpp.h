#ifndef DENSEGAUSS_E_PP_H
#define DENSEGAUSS_E_PP_H

#include "GaussianElimination.h"
class DenseGaussEPP : public GaussE
{
private:
    /* data */
public:
    DenseGaussEPP(/* args */);
    ~DenseGaussEPP();
    std::vector<double> solve(Matrix& A, std::vector<double>& b);
};

#endif
