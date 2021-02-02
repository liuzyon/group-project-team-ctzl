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

int main()
{
    cout << endl;
    cout << "--------------------------------------------" << endl;
    cout << "   Dense Matrix" << endl;
    cout << "--------------------------------------------" << endl;
//    cout << endl;
    // create A, b
    int rows = 3;
    int cols = 3;
    double input[9] = { 2, -3, 1, -2, 1, -2, 4, 0, 9};
    Matrix<double> A(rows, cols, input);
    double *b = new double [3]{ 2., 1., 5.};
    double *x = new double [3];
    cout << "--------------------------------" << endl;
    cout << "   Create Example A and b" << endl;
    cout << "--------------------------------" << endl;
    cout << "A created:" << endl;
    A.printMatrix();
    cout << endl;
    cout << "b created:" << endl;
    printVector(b, A.cols);
    cout << endl;

    Test<double> test;
    Solver<double> sv;



    cout << "<-------------------------------->" << endl;
    cout << "   DenseLUFactorisation Solve" << endl;
    cout << "<-------------------------------->" << endl;
    sv.DenseLUFactorisationSolve(A, b, x);

    cout << "x solved:" << endl;
    printVector(x, A.cols);

    test.testDense(A, b, x);


    delete[] b;
    delete[] x;
}

template <class T>
void printVector(T* vec, int size) {
    for (int i = 0; i < size; ++i)
    {
        cout << vec[i] << " ";
    }
    cout << endl;
}