//
// Created by Zhiyong Liu on 2021/1/26.
//

#include "DenseJacobi.h"
DenseJacobi::DenseJacobi(/* args */)
{
}

DenseJacobi::~DenseJacobi()
{
}

std::vector<double> DenseJacobi::solve(Matrix& A, std::vector<double>& b)
{
    // the output result of solver
    std::vector<double> x;
    const int n_iter = 10;
    int n = 0;
    int size = b.size();

    // initial guess of x is b
    x.reserve(size);
    for (int i = 0; i < size; i++)
    {
        x.push_back(b[i]);
    }

    std::vector<double> x_new;
    x_new.resize(size);

    while(n < n_iter)
    {

        for (int i = 0; i < size; i++)
        {
            int sum = 0;
            for (int j = 0; j < size; j++)
            {
                if (i != j)
                {
                    sum += A.values[i * size + j] * x[j];
                }
            }
            x_new[i] = (b[i] - sum) / A.values[i * size + i];
        }
        for (int i = 0; i < size; i++)
            {
                x[i] = x_new[i];
            }
        n++;
    }
    return x;
}