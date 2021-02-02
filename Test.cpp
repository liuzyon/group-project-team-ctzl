//
// Created by Zhiyong Liu on 2021/2/2.
//
#include <iostream>
#include "Test.h"
#include "Solver.h"
using namespace std;

template <class T>
Test<T>::Test()
{}

template <class T>
Test<T>::~Test()
{

}

template <class T>
bool Test<T>::testDense(Matrix<T>& A, T* b, T* x)
{
    cout << endl << "Check if Ax == b:" << endl;
    int size = A.cols;
    double *check_b = new double [A.cols];
    A.matVecMult(x, check_b);

    bool testPass = true;
    for (int i = 0; i < size; ++i)
    {
        if (check_b[i] != b[i]) testPass = false;
    }

    if (testPass) {
        cout << "Test Passed" << endl;
    } else {
        cout << "Test Failed" << endl;
    }

    return testPass;
}
