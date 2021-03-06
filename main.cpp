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

    // symmetirx matrix for Cholesky solver
    double* input1 = new double[rows * cols];
    for (int i = 0; i < rows * cols; i++)
    {
        input1[i] = rand() % 10 + 1;
    }
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            input1[i * cols + j] = input1[j * cols + i];
        }
    }
    // make diagnaolly dominant 
    for (int i = 0; i < rows; i++)
    {
        input1[i * (cols + 1)] = (rand() % 10 + 10) * 10;
    }
    Matrix<double> A1(rows, cols, input1);
    double *b = new double [rows];
    for (int i = 0; i < rows; i++) {
        b[i] = rand() % 5 + 1;
    }
    double *x = new double [rows];
    cout << "A created:" << endl;
    A.printMatrix();
    cout << endl;
    cout << "A1 created:" << endl;
    A1.printMatrix();
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

    printStartTag("Dense Cholesky Solve");
    sv.DenseCholeskySolve(A1, b, x);
    cout << "x solved:" << endl;
    printVector(x, A1.cols);
    test.test_result(A1, b, x);
    printEndTag();

    printStartTag("Dense GMRES Solve (only part of this algorithm finished) ");
    //// only part of DenseGMRES implemented
    // sv.DenseGMRES(A, b, x);
    printEndTag();

    delete[] b;
    delete[] x;
    delete[] input;
    delete[] input1;

    cout << endl;
    cout << "------------------------------------------------------" << endl;
    cout << "                   Sparse Matrix" << endl;
    cout << "------------------------------------------------------" << endl;
    cout << "--------------------------------" << endl;
    cout << "   Create Example A and b:" << endl;
    cout << "--------------------------------" << endl;
    // create A, b, x 
    int rows_sparse1 = 11;
    int cols_sparse1 = 11;
    int nnzs1 = 20;
    double values1[20] = { 100 ,100, 3, 100, 8, 7, 100,100,6,6,100,8,100,100,12,12,100,100,3,100  };
    int row_position1[12] = { 0, 1, 3, 5, 7 ,9,11,13,15,17,18,20};
    int col_index1[20] = { 0, 1, 10, 2, 6, 2, 3,4,5,4,5,2,6,7,8,7,8 ,9,1,10 };
    CSRMatrix<double> A1_sparse(rows_sparse1, cols_sparse1, nnzs1, values1, row_position1, col_index1);
    double* b1_sparse = new double[11];
    for (int i = 0; i < 11; i++) {
        b1_sparse[i] = i + 1;
    }

    double* x1_sparse = new double[11];
    cout << "A_sparse created:" << endl;
    A1_sparse.printMatrix();
    cout << endl;
    cout << "b_sparse created:" << endl;
    printVector(b1_sparse, A1_sparse.cols);

    Test<double> test_sparse;
    Solver<double> sv_sparse;
    cout << "--------------------------------" << endl;
    cout << "   Initialization completed!" << endl;
    cout << "--------------------------------" << endl;

    printStartTag("Sparse Jacobi Solve (using default tolerance) ");
    // can add another parameter as tolerance
    sv_sparse.sparse_jacobi_solver(A1_sparse, b1_sparse, x1_sparse);
    cout << "x_sparse solved:" << endl;
    printVector(x1_sparse, A1_sparse.cols);
    test_sparse.test_result(A1_sparse, b1_sparse, x1_sparse);
    printEndTag();

    printStartTag("Sparse Gauss-Seidel Solve (using default tolerance) ");
    // can add another parameter as tolerance
    sv_sparse.sparse_gauss_seidel_solver(A1_sparse, b1_sparse, x1_sparse);
    cout << "x_sparse solved:" << endl;
    printVector(x1_sparse, A1_sparse.cols);
    test_sparse.test_result(A1_sparse, b1_sparse, x1_sparse);
    printEndTag();


    printStartTag("Sparse Multigrid Solve (only work with odd rows/cols) ");
    sv_sparse.sparse_multigrid_solver(A1_sparse, b1_sparse, x1_sparse);
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

    delete[] b1_sparse;
    delete[] x1_sparse;
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