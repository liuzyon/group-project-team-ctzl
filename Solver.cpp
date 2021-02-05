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
void Solver<T>::DenseGaussESolve(Matrix<T> A, T b[], T* x)
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
void Solver<T>::DenseGaussEPPSolve(Matrix<T> A, T b[], T* x)
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

    T *L = new T[A.rows*A.cols]();
    T *U = new T[A.rows*A.cols]();

    // copy A to U
    for (int i = 0; i < A.rows*A.cols; ++i)
    {
        U[i] = A.values[i];
    }


    // Loop over each pivot row - don't need to consider the final row as a pivot
    for (int k = 0; k < A.rows-1; ++k)
    {
        // Loop over each equation below the pivot row - now we do need to consider the last row
        for (int i = k+1; i <A.rows ; ++i)
        {
            // Define the scaling factor outside the innermost
            // loop otherwise its value gets changed.
            T s = U[i*A.cols+k] / U[k*A.cols+k];
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
    std::vector<T> y;
    y.resize(size);
    for (int k = 0; k < size; ++k)
    {
        T s = 0;
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
        T s = 0;
        for (int j = k+1; j < size; ++j)
        {
            s = s + U[k*A.cols+j] * x[j];
        }
        x[k] = (y[k] - s) / U[k*A.cols+k];
    }

    delete[] L;
    delete[] U;
}

template <class T>
void Solver<T>::dense_multigrid_solver(Matrix<T>& A, T* b, T* x)
{
    
    if (A.rows % 2 == 0)
    {
        std::cerr << "            Caution!!!" << std::endl;
        std::cerr << "This case has not been implemented." << std::endl;
        std::cerr << "Ignored if passed in test." << std::endl;
    }
    // for the multigrid method, we need a restriction matrix to change from fine grid to coarse grid
    // and also a interpolation matrix to change from coarse grid back to the fine grid

    // the following code only works for input matrix with odd number of rows/cols
    // since I will define a coarse grid with double the size of patition in fine grid
    // meaning a 7 element vector will be changed to 3 etc.
    int fine_size = A.cols;
    int coarse_size = (fine_size - 1) / 2;

    // create the interpolation matrix fisrt
    // the interplation matrix has a fine_size as number of rows and coarse_size as number of cols

    // store the interpolation matrix as each col contain 1/2, 1, 1/2
    // see printMatrix()
    int elements = fine_size * coarse_size;
    double* values_inter = new double[elements];

    // fill with 0
    for(int i = 0; i < elements; i++)
    {
        values_inter[i] = 0;
    }

    // add in 1/2 1 at proper location
    values_inter[0] = 0.5;
    values_inter[elements - 1] = 0.5;
    // ignore first and last row, defined above
    for(int i = 1; i < fine_size - 1; i++)
    {
        if (i % 2 == 1)
            values_inter[(i - 1) / 2 + i * coarse_size] = 1;
        if (i % 2 == 0)
        {
            values_inter[(i / 2) + (i * coarse_size)] = 0.5;
            values_inter[(i / 2) - 1 + (i * coarse_size)] = 0.5;    
        }
    }

    // create the Interpolation matrix
    // the interplation matrix has a fine_size as number of rows and coarse_size as number of cols
    Matrix<double>* Inter = new Matrix<double>(fine_size, coarse_size, values_inter);

    // the Restrition matrix can be found as R = 1/2 * transpose of (I)
    // rows and cols will be swapped
    // fill with 0 first
    double* values_restr = new double[elements];
    for (int i = 0; i < elements; i++)
    {
        values_restr[i] = 0;
    }
    // add in 1/2 1/4 at proper location
    for(int i = 0; i < coarse_size; i++)
    {
        values_restr[i * (fine_size + 2)] = 0.25;
        values_restr[1 + i * (fine_size + 2)] = 0.5;
        values_restr[2 + i * (fine_size + 2)] = 0.25;
    }

    // create the restriction matrix 
    // the restriction matrix has coarse_size as number of rows and fine_size as number of cols
    Matrix<double>* Restr= new Matrix<double>(coarse_size, fine_size, values_restr);

    std::cout << "Interpolation martix used: ";
    Inter->printMatrix();
    std::cout << std::endl;
    std::cout << "Restriction martix used: ";
    Restr->printMatrix();
    std::cout << std::endl;

    // Setup complete

    double tol = 1;
    const int n_iter = 10;
    int n = 0;

    // set x to be same as b at the beginning
    for (int i = 0; i < fine_size; i++)
    {
        x[i] = b[i];
    }

    // stop when max iteration reached (now is 10), or a tolerance has been reached, now set to 1e-6
    while(tol > 1e-6 && n < n_iter)
    {
        // Now I will use n iterations in dense gauss seidel to flatten the error, here I use 2
        int n_pre_smooth = 0;
        const int n_iter_pre_smooth = 2;
        
        while (n_pre_smooth < n_iter_pre_smooth)
        {
            for (int i = 0; i < fine_size; i++)
            {
                double sum = 0;
                for (int j = 0; j < fine_size; j++)
                {
                    if (i != j)
                    {
                        sum += A.values[i * fine_size + j] * x[j];
                    }
                }
                // update x[i] after iterating each row
                x[i] = (b[i] - sum) * (1 / A.values[i * fine_size + i]);
            }
            n_pre_smooth++;
        }
        
        // std::cout << std::endl;
        // std::cout << "x after pre-smoothing: ";
        // for (int i = 0; i < fine_size; i++)
        // {
        //     std::cout << x[i] << "  ";
        // }

        double* result = new double[fine_size];
        double* error_fine = new double[fine_size];

        // compute the residual error after cycles of gauss-seidel
        A.matVecMult(x, result);
        for (int i = 0; i < fine_size; i++)
        {
            error_fine[i] = b[i] - result[i];
        }

        // Now use the restriction matrix to change the error in fine grid to error in coarse grid
        double* error_coarse = new double[coarse_size];
        Restr->matVecMult(error_fine, error_coarse);

        // Now change matrix A to coarse grid
        // To do this, I use R * A * I
        // matMatMult is a left multiplication
        Matrix<double>* temp_mat = new Matrix<double>(coarse_size, fine_size, false);
        Matrix<double>* A_coarse = new Matrix<double>(coarse_size, coarse_size, false);
        
        A.matMatMult(*Restr, *temp_mat);

        // temp_mat->printMatrix();

        Inter->matMatMult(*temp_mat, *A_coarse);

        // A_coarse->printMatrix();

        // using the matrix and error in coarse grid 
        // x_coarse to store the output
        double* x_coarse = new double[coarse_size];
        // calculate A_coarse * x_coarse = error_coarse
        // I use gaussian elimination solver to get the exact result
        // not using gauss-seidel as it might not converge
        DenseGaussESolve(*A_coarse, error_coarse, x_coarse);
        
        // use the interpolation matrix to change x_coarse back to fine grid and update x
        double* x_fine = new double[fine_size];
        Inter->matVecMult(x_coarse, x_fine);

        // update x
        for (int i = 0; i < fine_size; i++)
        {
            x[i] += x_fine[i];
        }
        
        // std::cout << std::endl;
        // std::cout << "x after restriction and interpolation: ";
        // for (int i = 0; i < fine_size; i++)
        // {
        //     std::cout << x[i] << "  ";
        // }     
        // std::cout << std::endl;

        // post smoothing using gauss seidel n times, here I use 2
        int n_post_smooth = 0;
        const int n_iter_post_smooth = 2;
        
        while (n_post_smooth < n_iter_post_smooth)
        {
            for (int i = 0; i < fine_size; i++)
            {
                double sum = 0;
                for (int j = 0; j < fine_size; j++)
                {
                    if (i != j)
                    {
                        sum += A.values[i * fine_size + j] * x[j];
                    }
                }
                // update x[i] after iterating each row
                x[i] = (b[i] - sum) * (1 / A.values[i * fine_size + i]);
            }
            n_post_smooth++;
        }

        // std::cout << "x after post-smmothing: ";
        // for (int i = 0; i < fine_size; i++)
        // {
        //     std::cout << x[i] << "  ";
        // }
        // std::cout << std::endl;

        // get the current residual error
        tol = 0;
        tol = get_error(A, b, x, fine_size, tol);
        // std::cout << "Current tolerance: " << tol << std::endl;
        n++;

        delete[] result;
        delete[] error_fine;
        delete[] error_coarse;
        delete[] x_coarse;
        delete[] x_fine;
        delete temp_mat;
        delete A_coarse;

    }   
    std::cout << std::endl;
    delete[] values_inter;
    delete[] values_restr;
    delete Inter;
    delete Restr;

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
void Solver<T>::sparse_multigrid_solver(CSRMatrix<T>& A, T* b, T* x)
{
    if (A.rows % 2 == 0)
    {
        std::cerr << "            Caution!!!" << std::endl;
        std::cerr << "This case has not been implemented." << std::endl;
        std::cerr << "Ignored if passed in test." << std::endl;
    }
    // for the multigrid method, we need a restriction matrix to change from fine grid to coarse grid
    // and also a interpolation matrix to change from coarse grid back to the fine grid

    // the following code only works for input matrix with odd number of rows/cols
    // since I will define a coarse grid with double the size of patition in fine grid
    int fine_size = A.cols;
    int coarse_size = (fine_size - 1) / 2;

    // create the interpolation matrix fisrt
    // the interplation matrix has a fine_size as number of rows and coarse_size as number of cols

    // store the sparse interpolation matrix as each col contain 1/2, 1, 1/2
    int nnzs = coarse_size * 3;
    double* values_inter = new double[nnzs];
    int* col_index_inter = new int[nnzs];
    int* row_position_inter = new int[fine_size + 1];

    int count = 0;
    for (int i = 0; i < nnzs; i++)
    {
        if(i % 3 == 0 || i % 3 == 2)
            values_inter[i] = 0.5;
        else 
            values_inter[i] = 1;

        if(i % 3 == 0 && i != 0)
            count++;
        // define col_index used in CSRMatrix construction
        col_index_inter[i] = count;
    }

    // define row position
    row_position_inter[0] = 0;
    row_position_inter[1] = 1;
    row_position_inter[fine_size] = nnzs;
    for (int i = 2; i < fine_size; i++)
    {
        if (i % 2 == 0)
            row_position_inter[i] = row_position_inter[i - 1] + 1;
        else 
            row_position_inter[i] = row_position_inter[i - 1] + 2;

    }

    // create the Interpolation matrix
    // the interplation matrix has a fine_size as number of rows and coarse_size as number of cols
    CSRMatrix<double>* Inter = new CSRMatrix<double>(fine_size, coarse_size, nnzs, values_inter, row_position_inter, col_index_inter);

    // the Restrition matrix can be found as R = 1/2 * transpose of (I)
    // So all values will be halved, nnzs will not change 
    // rows and cols will be swapped
    double* values_restr = new double[nnzs];
    int* col_index_restr = new int[nnzs];

    for (int i = 0; i < nnzs; i++)
    {
        values_restr[i] = values_inter[i] / 2;
    }

    col_index_restr[0] = 0;
    for (int i = 1; i < nnzs; i++)
    {
        if (i % 3 == 0)
            col_index_restr[i] = col_index_restr[i - 1];
        else   
            col_index_restr[i] = col_index_restr[i - 1] + 1;
    }


    int* row_position_restr = new int[coarse_size + 1];
    for (int i = 0; i < coarse_size + 1; i++)
    {
        row_position_restr[i] = i * 3;
    }

    // create the restriction matrix 
    // the restriction matrix has coarse_size as number of rows and fine_size as number of cols
    CSRMatrix<double>* Restr= new CSRMatrix<double>(coarse_size, fine_size, nnzs, values_restr, row_position_restr, col_index_restr);

    // Setup complete

    std::cout << "Interpolation martix used: ";
    Inter->printMatrix();
    std::cout << std::endl;
    std::cout << "Restriction martix used: ";
    Restr->printMatrix();
    std::cout << std::endl;

    std::cerr << "Solver unfinished yet." << std::endl << std::endl;

    // I stopped here as I dont have time to finish it
    // But with the interpolation and restriction matrices all set up in CSR format
    // the rest should be quite similar with that in the dense_multigrid solver

    // pre-smoothing using cycles of sparse_gauss_seidel_solver

    // get the residual error

    // using restriction matrix to change the error to a coaser grid error_coarse

    // change the matrix A to a coaser grid A_coarse bt R * A * I

    // calculate new x_coarse as A_coarse * x_coarse = error_coarse

    // using interpolation matrix to change x_coarse to x_fine

    // update x using x += x_fine;

    // post-smmothing using cycles of gauss-seidel

    // repeat all above processes if tolerance not reached
    
    delete[] values_restr;
    delete[] col_index_restr;
    delete[] row_position_restr;
    delete[] values_inter;
    delete[] col_index_inter;
    delete[] row_position_inter;
    
    delete Inter;
    delete Restr;
    
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




























