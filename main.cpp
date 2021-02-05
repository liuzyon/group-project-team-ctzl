#include <iostream>
#include <vector>
#include "Matrix.h"
#include "Matrix.cpp"
#include "CSRMatrix.h"
#include "CSRMatrix.cpp"
#include "Solver.h"
#include "Solver.cpp"
#include "Test.h"
#include "Test.cpp"

using namespace std;

template <class T>
void printVector(T* vector, int size);
void printStartTag(string s);
void printEndTag();

// all tests here assume a 1e-6 default tolerance to test whether calculated results are accurate enough
// all solvers which could take a user tolerance are now also using a 1e-6 default tolerance in order to pass the test
int main()
{
    cout << endl;
    cout << "------------------------------------------------------" << endl;
    cout << "                   Dense Matrix" << endl;
    cout << "------------------------------------------------------" << endl;
    cout << "--------------------------------" << endl;
    cout << "   Create Example A and b:" << endl;
    cout << "--------------------------------" << endl;
    // create A, b, x
    int rows = 11;
    int cols = 11;
    double* input = new double[rows * cols];
    for (int i = 0 ; i < rows * cols; i++)
    {
        input[i] = rand()%10 + 1;
    }
    // make diagnaolly dominant 
    for (int i = 0 ; i < rows; i++)
    {
        input[i * (cols + 1)] = (rand()%10 + 10) * 10;
    }
    Matrix<double> A(rows, cols, input);
    double *b = new double [rows]{1, 4, 6, 3, 5, 2, 6, 8, 6, 3, 9};
    double *x = new double [rows];
    cout << "A created:" << endl;
    A.printMatrix();
    cout << endl;
    cout << "b created:" << endl;
    printVector(b, A.cols);

    Test<double> test;
    Solver<double> sv;
    cout << "--------------------------------" << endl;
    cout << "   Initialization completed!" << endl;
    cout << "--------------------------------" << endl;

    printStartTag("Dense Jacobi Solve (using default tolerance) ");
    // can add another parameter as tolerance
    sv.dense_jacobi_solver(A, b, x);
    cout << "x solved:" << endl;
    printVector(x, A.cols);
    test.test_result(A, b, x);
    printEndTag();

    printStartTag("Dense Gauss-Seidel Solve (using default tolerance) ");
    // can add another parameter as tolerance
    sv.dense_gauss_seidel_solver(A, b, x);
    cout << "x solved:" << endl;
    printVector(x, A.cols);
    test.test_result(A, b, x);
    printEndTag();

    printStartTag("Dense Gaussian Elimination Solve");
    sv.DenseGaussESolve(A, b, x);
    cout << "x solved:" << endl;
    printVector(x, A.cols);
    test.test_result(A, b, x);
    printEndTag();

    printStartTag("Dense Gaussian Elimination PP Solve");
    sv.DenseGaussEPPSolve(A, b, x);
    cout << "x solved:" << endl;
    printVector(x, A.cols);
    test.test_result(A, b, x);
    printEndTag();

    printStartTag("Dense LU Factorisation Solve");
    sv.DenseLUFactorisationSolve(A, b, x);
    cout << "x solved:" << endl;
    printVector(x, A.cols);
    test.test_result(A, b, x);
    printEndTag();

    printStartTag("Dense Multigrid Solve (only work with odd rows/cols)");
    sv.dense_multigrid_solver(A, b, x);
    cout << "x solved:" << endl;
    printVector(x, A.cols);
    test.test_result(A, b, x);
    printEndTag();

    printStartTag("Dense GMRES Solve (only part of this algorithm finished) ");
    // only part of DenseGMRES implemented
    sv.DenseGMRES(A, b, x);
    printEndTag();

    delete[] b;
    delete[] x;
    delete[] input;


    cout << endl;
    cout << "------------------------------------------------------" << endl;
    cout << "                   Sparse Matrix" << endl;
    cout << "------------------------------------------------------" << endl;
    cout << "--------------------------------" << endl;
    cout << "   Create Example A and b:" << endl;
    cout << "--------------------------------" << endl;
    // create A, b, x
    int rows_sparse = 5;
    int cols_sparse = 5;
    int nnzs = 12;
    double values[12] = {10 ,2 ,3, 40, 5, 6, 70, 8, 9, 5, 12, 8};
    int row_position[6] = {0, 2, 5, 9, 11, 12};
    int col_index[12] = {0, 3, 0, 1, 3, 0, 2, 3, 4, 2, 3, 4};
    CSRMatrix<double> A_sparse(rows_sparse, cols_sparse, nnzs, values, row_position, col_index);
    double *b_sparse = new double [5]{5, 3, 8, 10, 4};
    double *x_sparse = new double [5];
    cout << "A_sparse created:" << endl;
    A_sparse.printMatrix();
    cout << endl;
    cout << "b_sparse created:" << endl;
    printVector(b_sparse, A_sparse.cols);

    Test<double> test_sparse;
    Solver<double> sv_sparse;
    cout << "--------------------------------" << endl;
    cout << "   Initialization completed!" << endl;
    cout << "--------------------------------" << endl;

    printStartTag("Sparse Jacobi Solve (using default tolerance) ");
    // can add another parameter as tolerance
    sv_sparse.sparse_jacobi_solver(A_sparse, b_sparse, x_sparse);
    cout << "x_sparse solved:" << endl;
    printVector(x_sparse, A_sparse.cols);
    test_sparse.test_result(A_sparse, b_sparse, x_sparse);
    printEndTag();

    printStartTag("Sparse Gauss-Seidel Solve (using default tolerance) ");
    // can add another parameter as tolerance
    sv_sparse.sparse_gauss_seidel_solver(A_sparse, b_sparse, x_sparse);
    cout << "x_sparse solved:" << endl;
    printVector(x_sparse, A_sparse.cols);
    test_sparse.test_result(A_sparse, b_sparse, x_sparse);
    printEndTag();

    printStartTag("Sparse Multigrid Solve (only work with odd rows/cols) ");
    sv_sparse.sparse_multigrid_solver(A_sparse, b_sparse, x_sparse);
    printEndTag();



    cout << endl;
    cout << "------------------------------------------------------" << endl;
    cout << "                   Test CSRMatrix multiply" << endl;
    cout << "------------------------------------------------------" << endl;

    int rows_a = 3;
    int cols_a = 3;
    int nnzs_a = 6;
    double a_values[12] = {1, 2, 3, 4, 5, 6};
    int a_row_position[6] = {0, 2, 3, 6};
    int a_col_index[12] = {0, 2, 2, 0, 1, 2};
    CSRMatrix<double> a_sparse(rows_a, cols_a, nnzs_a, a_values, a_row_position, a_col_index);
    std::cout << "Right Matrix:" << endl;
    a_sparse.printMatrix();

    std::cout << std::endl;

    int rows_c = 3;
    int cols_c = 3;
    int nnzs_c = 6;
    double c_values[12] = {1, 2, 3, 4, 5, 6};
    int c_row_position[6] = {0, 2, 3, 6};
    int c_col_index[12] = {0, 2, 2, 0, 1, 2};
    CSRMatrix<double> c_sparse(rows_c, cols_c, nnzs_c, c_values, c_row_position, c_col_index);
    std::cout << "Left Matrix:" << endl;
    c_sparse.printMatrix();

    printStartTag("Call the matMatMult for two sparse matrices");
    CSRMatrix<double> output(3, 3, -1, false);
    a_sparse.matMatMult(c_sparse, output);

    std::cout << "Output Matrix:" << endl;
    output.printMatrix();
    printEndTag();

    delete[] b_sparse;
    delete[] x_sparse;
}



template <class T>
void printVector(T* vec, int size) {
    cout << "[";
    for (int i = 0; i < size; ++i)
    {
        cout << vec[i];
        if (i != size-1) cout << ", ";
    }
    cout << "]" << endl;
}

void printStartTag(string s) {
    cout << endl << endl;
    cout << "<------------------------------------------" << endl;
    cout << "\t" << s << endl;
    cout << "<------------------------------------------" << endl;
}

void printEndTag() {
    cout << "---------------------->" << endl;
    cout << "   End" << endl;
    cout << "---------------------->" << endl;
}