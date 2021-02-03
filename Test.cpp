//
// Created by Zhiyong Liu on 2021/2/2.
//
#include <iostream>
#include "Test.h"
#include "Solver.h"
#include <math.h>
using namespace std;

template <class T>
Test<T>::Test()
{}

template <class T>
Test<T>::~Test()
{}

template <class T>
// test has a default tolerance of 1e-6
bool Test<T>::testDense(Matrix<T>& A, T* b, T* x, double tol)
{
    cout << endl << "Check if Ax == b:" << endl;
    int size = A.cols;
    double *check_b = new double [A.cols];
    A.matVecMult(x, check_b);

    bool testPass = false;
    double error = 0;
    for (int i = 0; i < size; ++i)
    {
        error += pow(check_b[i] - b[i], 2);
    }

    if (sqrt(error) < tol)
        testPass = true;

    if (testPass) {
        cout << "Test Passed" << endl;
    } else {
        cout << "Test Failed" << endl;
    }

    delete[] check_b;

    return testPass;
}
