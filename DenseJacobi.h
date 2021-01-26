//
// Created by Zhiyong Liu on 2021/1/26.
//

#ifndef DENSEJACOBI_H
#define DENSEJACOBI_H

#include "Jacobi.h"
class DenseJacobi : public Jacobi
{
private:
    /* data */
public:
    DenseJacobi(/* args */);
    ~DenseJacobi();
    std::vector<double> solve(Matrix& A, std::vector<double>& b);
};

#endif


