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