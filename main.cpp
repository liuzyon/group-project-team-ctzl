#include <iostream>
#include <vector>
#include "Matrix.h"
#include "Matrix.cpp"
#include "CSRMatrix.h"
#include "Solver.h"
#include "Solver.cpp"
#include "Test.h"
#include "Test.cpp"

using namespace std;

template <class T>
void printVector(T* vector, int size);
void printStartTag(string s);
void printEndTag();

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
    int rows = 3;
    int cols = 3;
    double input[9] = { 2, -3, 1, -2, 1, -2, 4, 0, 9};
    Matrix<double> A(rows, cols, input);
    double *b = new double [3]{ 2., 1., 5.};
    double *x = new double [3];
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

    printStartTag("Dense Gaussian Elimination Solve");
    sv.DenseGaussESolve(A, b, x);
    cout << "x solved:" << endl;
    printVector(x, A.cols);
    test.testDense(A, b, x);
    printEndTag();

    printStartTag("Dense Gaussian Elimination PP Solve");
    sv.DenseGaussEPPSolve(A, b, x);
    cout << "x solved:" << endl;
    printVector(x, A.cols);
    test.testDense(A, b, x);
    printEndTag();

    printStartTag("Dense LU Factorisation Solve");
    sv.DenseLUFactorisationSolve(A, b, x);
    cout << "x solved:" << endl;
    printVector(x, A.cols);
    test.testDense(A, b, x);
    printEndTag();


    delete[] b;
    delete[] x;

    cout << endl;
    cout << "------------------------------------------------------" << endl;
    cout << "          Diaganolly Dominant Dense Matrix" << endl;
    cout << "------------------------------------------------------" << endl;
    cout << "--------------------------------" << endl;
    cout << "   Create Example A and b:" << endl;
    cout << "--------------------------------" << endl;
    // create A, b, x
    int rows_2 = 4;
    int cols_2 = 4;
    double input_2[16] = {10, 2, 3, -5, 1, 14, 3, 2, -1, 4, 16, -4, 5, 4, 3, 21};
    Matrix<double> A_2(rows_2, cols_2, input_2);
    double *b_2 = new double [4]{3, 2, 5, 4};
    double *x_2 = new double [4];
    cout << "A_2 created:" << endl;
    A_2.printMatrix();
    cout << endl;
    cout << "b_2 created:" << endl;
    printVector(b_2, A_2.cols);

    Test<double> test_2;
    Solver<double> sv_2;
    cout << "--------------------------------" << endl;
    cout << "   Initialization completed!" << endl;
    cout << "--------------------------------" << endl;

    printStartTag("Dense Jacobi Solve (using default tolerance) ");
    sv_2.dense_jacobi_solver(A_2, b_2, x_2, 1e-6);
    cout << "x_2 solved:" << endl;
    printVector(x_2, A_2.cols);
    test_2.testDense(A_2, b_2, x_2);
    printEndTag();

    printStartTag("Dense Gauss-Seidel Solve (using default tolerance) ");
    sv_2.dense_gauss_seidel_solver(A_2, b_2, x_2);
    cout << "x_2 solved:" << endl;
    printVector(x_2, A_2.cols);
    test_2.testDense(A_2, b_2, x_2);
    printEndTag();

    delete[] b_2;
    delete[] x_2;
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