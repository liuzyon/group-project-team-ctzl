#include "Solver.h"
#include "Matrix.h"
#include "CSRMatrix.h"
#include <math.h>
#include <iostream>

template <class T>
Solver<T>::Solver(/* args */)
{
}

template <class T>
Solver<T>::~Solver()
{
}

// this is a private function used in dense/sparse Jacobi and Gauss-Seidel
// to calculate the error after 1 iteration and then compare with user tolerance
template <class T>
double Solver<T>::get_error(Matrix<T>& A, T* b, T* x, int size, double error)
{
    // calculate error between expected x and our A*x
    // create a output pointer to store the result
    T* output = new T[size];
    // matvecMult is a vitual funtion implemented in both dense and sparse matrix
    // depends on different matrix as argument, different versions will be called
    A.matVecMult(x, output);
    // compute the error
    for (int i = 0; i < size; i++)
    {
        // calculate the norm
        error += pow((output[i] - b[i]),2);
    }
    delete[] output;

    // calculate the norm
    return sqrt(error);
}

//Solver starts from below
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

template <class T>
// there is a default user tolerance, can be changed by user
void Solver<T>::dense_jacobi_solver(Matrix<T>& A, T* b, T* x, double tol)
{
    // max iteration if tolerance not reached
    const int n_iter = 1000;
    int n = 0;
    int size = A.cols;

    // initial guess of x is b
    for (int i = 0; i < size; i++)
    {
        x[i] = b[i];
    }

    while(n < n_iter)
    {
        // x new here stores and updates all x values after 1 iteration
        T* x_new = new T[size];
        // error keeps track of the error after 1 iteration
        double error = 0;

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
        // update x using x_new
        for (int i = 0; i < size; i++)
        {
            x[i] = x_new[i];
        }
        delete[] x_new;

        // get the current error
        error = get_error(A, b, x, size, error);

        // if below the tolerance, ends the iteration
        if (error < tol)
            return;
        // else iterates
        else
            n++;
    }
}

template <class T>
void Solver<T>::dense_gauss_seidel_solver(Matrix<T>& A, T* b, T* x, double tol)
{
    // max iteration if tolerance not reached
    const int n_iter = 1000;
    int n = 0;
    int size = A.cols;

    // initial guess of x is b
    for (int i = 0; i < size; i++)
    {
        x[i] = b[i];
    }

    while(n < n_iter)
    {
        double error = 0;

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
            // update x[i] after iterating each row
            x[i] = (b[i] - sum) * (1 / A.values[i * size + i]);
        }

        // get the current error
        error = get_error(A, b, x, size, error);

        // if below the tolerance, ends the iteration
        if (error < tol)
            return;
        // else iterates
        else
            n++;
    }
}


template <class T>
void Solver<T>::DenseGaussESolve(Matrix<T> &A, T* b, T* x)
{
    // This method is direct calculate the answer of Ax=b.
    // This solver is only for the matix which exist the only answer x and no 0 on diagnal entries.
    // The computaional cost of Gaussian Elimination is O(n^3).
    //and we will also highlight issues that arise when the value is non-zero but small in pivoting.
    //make sure the sizes of matixs are what we can calculate
    int size = A.cols;

    if (A.cols != A.rows) {
        std::cerr << "Input dimensions for matrices don't match" << std::endl;
        return;
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
            for (int j = k; j < size; j++) {
                A.values[i * A.cols + j] = A.values[i * A.cols + j] - s * A.values[k * A.cols + j];
            }
            // and update the corresponding entry of b
            b[i] = b[i] - s * b[k];
        }
    }

    //backsubstitude

    //start at the end (row n-1) and work backwards
    for (int k = size - 1; k > -1; k--) {
        s = 0;
        for (int j = k + 1; j < size; j++) {
            //sum of Akj*xj from k+1 to n
            s = s + A.values[k * A.cols + j] * x[j];
        }
        //xk = (1/akk)*(bk - sum of Akj*xj from k+1 to n)
        x[k] = (b[k] - s) / A.values[k * A.cols + k];
    }

}

template <class T>
void Solver<T>::DenseGaussEPPSolve(Matrix<T> &A, T* b, T* x)
{
    // This method is direct calculate the answer of Ax=b.
    // This solver is only for the matix which exist the only answer x.
    // The computaional cost of Gaussian Elimination is O(n^3).
    //and we fix the issue we meet in DenseGaussE class we use partial pivoting.
    //This is an upgrade of Gaussian Elimination method.
    int size = A.cols;
    //make sure the sizes of matixs are what we can calculate
    if (A.cols != A.rows) {
        std::cerr << "Input dimensions for matrices don't match" << std::endl;
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
            for (int j = k; j < size; j++) {
                A.values[i * A.cols + j] = A.values[i * A.cols + j] - s * A.values[k * A.cols + j];
            }
            // and update the corresponding entry of b
            b[i] = b[i] - s * b[k];
        }
    }

    //backsubstitude

    //start at the end (row n-1) and work backwards
    for (int k = size - 1; k > -1; k--) {
        s = 0;
        for (int j = k + 1; j < size; j++) {
            //sum of Akj*xj from k+1 to n
            s = s + A.values[k * A.cols + j] * x[j];
        }
        //xk = (1/akk)*(bk - sum of Akj*xj from k+1 to n)
        x[k] = (b[k] - s) / A.values[k * A.cols + k];
    }

}

 
template <class T>
void Solver<T>::DenseLUFactorisationSolve(Matrix<T>& A, T* b, T* x)
{
    int size = A.cols;
    //make sure the sizes of matixs are what we can calculate
    if (A.cols != A.rows) {
        std::cerr << "Input dimensions for matrices don't match" << std::endl;
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
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// solver for sparse matrix
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

template <class T>
void Solver<T>::sparse_jacobi_solver(CSRMatrix<T>& A, T* b, T* x, double tol)
{
    // max iteration if tolerance not reached
    const int n_iter = 1000;
    int n = 0;
    int size = A.cols;

    // initial guess of x is b
    for (int i = 0; i < size; i++)
    {
        x[i] = b[i];
    }

    while(n < n_iter)
    {
        // x new here stores and updates all x values after 1 iteration
        T* x_new = new T[size];
        // error keeps track of the error after 1 iteration
        double error = 0;

        // loop over row
        for (int i = 0; i < size; i++)
        {
            double sum = 0;
            double dia_val = 0;
            // loop over all elements in row i
            for (int val_index = A.row_position[i]; val_index < A.row_position[i+1]; val_index++)
            {
                // get the column index of that element
                int c_index = A.col_index[val_index];
                // check if it is in the diaganol
                if (i != c_index)
                {
                    sum += A.values[val_index] * x[c_index];
                }
                // get the diagnaol value
                else
                {
                    dia_val = A.values[val_index];
                }
            }
            x_new[i] = (b[i] - sum) * (1 / dia_val);
        }
        for (int i = 0; i < size; i++)
            {
                x[i] = x_new[i];
            }
        delete[] x_new;

        // get the current error
        error = get_error(A, b, x, size, error);

        // if below the tolerance, ends the iteration
        if (error < tol)
            return;
        // else iterates
        else
            n++;
    }
}



template <class T>
void Solver<T>::sparse_gauss_seidel_solver(CSRMatrix<T>& A, T* b, T* x, double tol)
{
    // max iteration if tolerance not reached
    const int n_iter = 1000;
    int n = 0;
    int size = A.cols;

    // initial guess of x is b
    for (int i = 0; i < size; i++)
    {
        x[i] = b[i];
    }

    while(n < n_iter)
    {
        double error = 0;

        // loop over row
        for (int i = 0; i < size; i++)
        {
            double sum = 0;
            double dia_val = 0;
            // loop over all elements in row i
            for (int val_index = A.row_position[i]; val_index < A.row_position[i+1]; val_index++)
            {
                // get the column index of that element
                int c_index = A.col_index[val_index];
                // check if it is in the diaganol
                if (i != c_index)
                {
                    sum += A.values[val_index] * x[c_index];
                }
                // get the diagnaol value
                else
                {
                    dia_val = A.values[val_index];
                }
            }
            // update before getting into another row
            x[i] = (b[i] - sum) * (1 / dia_val);
        }

        // get the current error
        error = get_error(A, b, x, size, error);

        // if below the tolerance, ends the iteration
        if (error < tol)
            return;
        // else iterates
        else
            n++;
    }
}

template <class T>
void Solver<T>::DenseGMRES(Matrix<T> &A, T* b, T* x)
{

    //make sure the sizes of matixs are what we can calculate
    if (A.cols != A.rows) {
        std::cerr << "Input dimensions for matrices don't match" << std::endl;
    }


    int n = A.cols; // dimension of each vector: n

    Matrix<T> A_copy(A);

    // Arnoldi algorithm:
    // 由向量组b, Ab, A^2b,..., A^k-1v生成的向量空间称为k维Krylov子空间，记作K(k)
    // 求一组与之等价的规范正交向量组
    // 矩阵A的各列用来存放输入的列向量，b为输入的列向量，k为迭代次数

    // 最大迭代步数IterMax
    int k = 4;

//    std::unique_ptr<T*[]> K(new T*[k]);     //K向量组
    T **K_N = new T*[k];
    for (int i = 0; i < k; ++i)
    {
        K_N[i] = new T[n];
    }
//    std::unique_ptr<T*[]> K_N(new T*[k]);   //与K等价的规范正交向量组
//    T *v0 = new T[n];
//    std::unique_ptr<T[]> v0(new T[n]);

    T **H = new T*[k];
    for (int i = 0; i < k; ++i)
    {
        H[i] = new T[k+1];
    }



    // 初值x0
    T *x0 = new T[n];
    // 将初值x0设为b0
    for (int i = 0; i < n; ++i)
    {
        x0[i] = b[0];
    }

    T *Ax0 = new T[n];
    A.matVecMult(x0, Ax0);
    T *r0 = new T[n];
    for (int i = 0; i < n; ++i)
    {
        r0[i] = b[i] - Ax0[0];
    }

    T beta = calMagnitude(r0, n);
    T *v0 = new T[n];
    unitization(r0, v0, n);
    // v1
    K_N[0] = v0;

    for (int j = 0; j < k; ++j)
    {
        // A_vj
        T *w = new T[n];
        A.matVecMult(K_N[j], w);


        for (int i = 0; i <= j; ++i)
        {
            // A_vj和每个已产生的vi作内积: hij

            H[i][j] = innerProduct(w, K_N[i], n); // 计算(Avj, vi)

            // 计算 hij * vi
            T* pro_vec = new T[n];
            for (int l = 0; l < n; ++l)
            {
                pro_vec[l] = K_N[i][l] * H[i][j];
            }

            // w = w - hij * vj
            for (int l = 0; l < n; ++l)
            {
                w[l] -= pro_vec[l];
            }

            delete[] pro_vec;
        }

        // hj+1,j
        T hnext = calMagnitude(w, n);
        if (hnext == 0) {
            // m = j, break;
            break;
        }

        unitization(w, K_N[j+1], n);

        // 检测是否停机
        T* rj = new T[n];
        // 如何计算残量rm的范数


        for (int i = 0; i < j-1; ++i)
        {

        }




        delete[] w;
    }


//    for (int i = 0; i < k; ++i)
//    {
//        for (int j = 0; j < n; ++j)
//        {
//            std::cout << K_N[i][j] << " ";
//        }
//        std::cout << std::endl;
//    }

//    for (int i=0; i<k; i++) {
//        delete[] K_N[i];
//    }
    delete[] K_N;
    delete[] v0;
}


template <class T>
T Solver<T>::calMagnitude(T *input, int size)
{
    T sum = 0;
    for (int i = 0; i < size; ++i)
    {
        sum += pow(input[i], 2);
    }
    T magnitude = pow(sum, 0.5);
}

template <class T>
void Solver<T>::unitization(T* input, T* output, int size)
{
    T sum = 0;
    for (int i = 0; i < size; ++i)
    {
        sum += pow(input[i], 2);
    }
    T magnitude = pow(sum, 0.5);

    for (int i = 0; i < size; ++i)
    {
        output[i] = input[i] / magnitude;
    }
}


template <class T>
T Solver<T>::innerProduct(T *one, T *two, int size)
{
    T sum = 0;
    for (int i = 0; i < size; ++i)
    {
        sum += one[i] * two[i];
    }
    return sum;
}




























