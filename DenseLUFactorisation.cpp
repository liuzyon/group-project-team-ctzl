//
// Created by Zhiyong Liu on 2021/1/28.
//
#include <iostream>
#include "DenseLUFactorisation.h"

DenseLUFactorisation::DenseLUFactorisation()
{}

DenseLUFactorisation::~DenseLUFactorisation()
{

}

std::vector<double> DenseLUFactorisation::solve(Matrix& A, std::vector<double>& b)
{
    std::vector<double> x;
    x.resize(b.size());
    //make sure the sizes of matixs are what we can calculate
    if (A.cols != A.rows || A.cols != b.size()) {
        std::cerr << "Input dimensions for matrices don't match" << std::endl;
        return x;
    }

    double *L = new double[A.rows*A.cols]();
    double *U = new double[A.rows*A.cols]();

    // copy A to U
    for (int i = 0; i < A.rows; ++i)
    {
        for (int j = 0; j < A.cols; ++j)
        {
            U[i*A.cols+j] = A.values[i*A.cols+j];
        }
    }


    // Loop over each pivot row - don't need to consider the final row as a pivot
    for (int k = 0; k < A.rows-1; ++k)
    {
        // Loop over each equation below the pivot row - now we do need to consider the last row
        for (int i = k+1; i <A.rows ; ++i)
        {
            // Define the scaling factor outside the innermost
            // loop otherwise its value gets changed.
            double s = U[i*A.cols+k] / U[k*A.cols+k];
            for (int j = k; j <A.cols ; ++j)
            {
                U[i*A.cols+j] = U[i*A.cols+j] - s*U[k*A.cols+j];
            }
            // store the scaling factors which make up the lower tri matrix
            L[i*A.cols+k] = s;
        }
    }
    // remember to add in the ones on the main diagonal to L
    for (int i = 0; i < A.rows; ++i)
    {
        for (int j = 0; j < A.cols; ++j)
        {
            if (i == j) {
                L[i*A.cols+j]++;
            }
        }
    }


    //Forward substitution
    // A = L, b = b, x = y
    std::vector<double> y;
    int size = b.size();
    y.resize(size);
    for (int k = 0; k < size; ++k)
    {
        double s = 0;
        for (int j = 0; j < k; ++j)
        {
            s = s + L[k*A.cols+j] * y[j];
        }
        y[k] = (b[k] - s) / L[k*A.cols+k];
    }

    //Backward substitution
    // A = U, b = y, x = x
    for (int k = size-1; k > -1 ; --k)
    {
        double s = 0;
        for (int j = k+1; j < size; ++j)
        {
            s = s + U[k*A.cols+j] * x[j];
        }
        x[k] = (y[k] - s) / U[k*A.cols+k];
    }

    delete[] L;
    delete[] U;
    return x;
}
