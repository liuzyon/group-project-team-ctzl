#include "DenseGaussSeidel.h"

DenseGaussSeidel::DenseGaussSeidel(/* args */)
{
}

DenseGaussSeidel::~DenseGaussSeidel()
{
}

std::vector<double> DenseGaussSeidel::solve(Matrix& A, std::vector<double>& b)
{
    // the output result of solver
    std::vector<double> x;
    const int n_iter = 100;
    int n = 0;
    int size = b.size();

    // initial guess of x is b
    x.reserve(size);
    for (int i = 0; i < size; i++)
    {
        x.push_back(b[i]);
    }

    while(n < n_iter)
    {

        for (int i = 0; i < size; i++)
        {
            double sum = 0;
            for (int j = 0; j < size; j++)
            {
                if (i != j)
                {
                    sum += A.values[i * size + j] * x[j];
                }
            }
            x[i] = (b[i] - sum) / A.values[i * size + i];
        }
        n++;
    }
    return x;
}