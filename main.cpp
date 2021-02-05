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
    int rows = 4;
    int cols = 4;
    double input[16] = {10, 2, 3, -5, 1, 14, 3, 2, -1, 4, 16, -4, 5, 4, 3, 21};
    Matrix<double> A(rows, cols, input);
    double *b = new double [4]{ 2., 1., 5., 6.};
    double *x = new double [4];
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

    delete[] b;
    delete[] x;


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

    // sv.DenseGMRES(A, b, x);

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