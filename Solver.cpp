
#include "Solver.h"
#include "Matrix.h"
#include <iostream>

Solver::Solver(/* args */)
{
}

Solver::~Solver()
{
}

std::vector<double> Solver::DenseGaussESolve(Matrix &A, std::vector<double> &b)
{
    // This method is direct calculate the answer of Ax=b.
    // This solver is only for the matix which exist the only answer x and no 0 on diagnal entries.
    // The computaional cost of Gaussian Elimination is O(n^3).
    //and we will also highlight issues that arise when the value is non-zero but small in pivoting.
    std::vector<double> x;
    x.resize(b.size());
    //make sure the sizes of matixs are what we can calculate
    if (A.cols != A.rows || A.cols != b.size()) {
        std::cerr << "Input dimensions for matrices don't match" << std::endl;
        return x;
    }
    //triangularise
    // covert A into upper triangluar form through row operations.
    //The same row operations are performed on the vector b.
    double s;
    // Loop over each pivot row all but the last row which we will never need to use as a pivot
    for (int k = 0; k < A.rows -1; k++) {
        //Loop over each row below the pivot row, including the last row which we do need to update
        for (int i = k + 1; i < A.rows; i++) {
            //Define the scaling factor for this row outside the innermost loop otherwise
            //There's also a performance saving from not recomputing things when not strictly necessary
            s = A.values[i * A.cols + k] / A.values[k * A.cols + k];
            //Update the current row of A by looping over the column j
            //start the loop from k as we can assume the entries before this are already zero
            for (int j = k; j < b.size(); j++) {
                A.values[i * A.cols + j] = A.values[i * A.cols + j] - s * A.values[k * A.cols + j];
            }
            // and update the corresponding entry of b
            b[i] = b[i] - s * b[k];
        }
    }

    //backsubstitude

    //start at the end (row n-1) and work backwards
    for (int k = b.size() - 1; k > -1; k--) {
        s = 0;
        for (int j = k + 1; j < b.size(); j++) {
            //sum of Akj*xj from k+1 to n
            s = s + A.values[k * A.cols + j] * x[j];
        }
        //xk = (1/akk)*(bk - sum of Akj*xj from k+1 to n)
        x[k] = (b[k] - s) / A.values[k * A.cols + k];
    }

    return x;
}

std::vector<double> Solver::DenseGaussEPPSolve(Matrix &A, std::vector<double> &b)
{
    // This method is direct calculate the answer of Ax=b.
    // This solver is only for the matix which exist the only answer x.
    // The computaional cost of Gaussian Elimination is O(n^3).
    //and we fix the issue we meet in DenseGaussE class we use partial pivoting.
    //This is an upgrade of Gaussian Elimination method.
    std::vector<double> x;
    x.resize(b.size());
    //make sure the sizes of matixs are what we can calculate
    if (A.cols != A.rows || A.cols != b.size()) {
        std::cerr << "Input dimensions for matrices don't match" << std::endl;
        return x;
    }
    //triangularise
    // covert A into upper triangluar form through row operations.
    //The same row operations are performed on the vector b.
    double s;
    int kmax;
    // Loop over each pivot row all but the last row which we will never need to use as a pivot
    for (int k = 0; k < A.rows - 1; k++) {
        //Swap rows so we are always dividing through by the largest magnitude number.
        //initiatise kmax with the current pivot row (k)
        kmax = k;
        //loop over all entries below the pivot and select the k with the largest abs value
        for (int i = k + 1; i < A.rows; i++) {
            if (abs(A.values[kmax * A.cols + k]) < abs(A.values[i * A.cols + k])) {
                kmax = i;
            }

        }
        //and swap the current pivot row (k) with the row with the largest abs value below the pivot
        std::vector<double> iA;
        iA.resize(A.cols);
        for (int j = 0; j < A.cols; j++) {
            iA[j] = A.values[kmax * A.cols + j];
        }
        double ib = b[kmax];

        for (int j = 0; j < A.cols; j++) {
            A.values[kmax * A.cols + j] = A.values[k * A.cols + j];
        }
        b[kmax] = b[k];
        for (int j = 0; j < A.cols; j++) {
            A.values[k * A.cols + j] = iA[j];
        }
        b[k] = ib;


        for (int i = k + 1; i < A.rows; i++) {


            //Define the scaling factor for this row outside the innermost loop otherwise
            //There's also a performance saving from not recomputing things when not strictly necessary
            s = A.values[i * A.cols + k] / A.values[k * A.cols + k];
            //Update the current row of A by looping over the column j
            //start the loop from k as we can assume the entries before this are already zero
            for (int j = k; j < b.size(); j++) {
                A.values[i * A.cols + j] = A.values[i * A.cols + j] - s * A.values[k * A.cols + j];
            }
            // and update the corresponding entry of b
            b[i] = b[i] - s * b[k];
        }
    }

    //backsubstitude

    //start at the end (row n-1) and work backwards
    for (int k = b.size() - 1; k > -1; k--) {
        s = 0;
        for (int j = k + 1; j < b.size(); j++) {
            //sum of Akj*xj from k+1 to n
            s = s + A.values[k * A.cols + j] * x[j];
        }
        //xk = (1/akk)*(bk - sum of Akj*xj from k+1 to n)
        x[k] = (b[k] - s) / A.values[k * A.cols + k];
    }

    return x;

}

std::vector<double> Solver::DenseGaussSeidelSolve(Matrix &A, std::vector<double> &b)
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

std::vector<double> Solver::DenseJacobiSolve(Matrix &A, std::vector<double> &b)
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

    std::vector<double> x_new;
    x_new.resize(size);

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

std::vector<double> Solver::DenseLUFactorisationSolve(Matrix &A, std::vector<double> &b)
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
