//
// Created by Zhiyong Liu on 2021/2/2.
//

#ifndef ASSESSMENT1_TEST_H
#define ASSESSMENT1_TEST_H

#include "Matrix.h"

template <class T>
class Test
{
public:
    Test();

    virtual ~Test();

    // test has a default tolerance of 1e-6
    bool testDense(Matrix<T>& A, T* b, T* x, double tol = 1e-6);
};
#endif //ASSESSMENT1_TEST_H
