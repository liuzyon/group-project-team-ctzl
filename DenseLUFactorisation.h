//
// Created by Zhiyong Liu on 2021/1/28.
//

#ifndef DENSELUFACTORISATION_H
#define DENSELUFACTORISATION_H
#include "LUFactorisation.h"
class DenseLUFactorisation : public LUFactorisation
{
public:
    DenseLUFactorisation();

    ~DenseLUFactorisation();
    std::vector<double> solve(Matrix& A, std::vector<double>& b);
};

#endif //DENSELUFACTORISATION_H

