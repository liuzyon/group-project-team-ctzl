//
// Created by Zhiyong Liu on 2021/1/27.
//

#ifndef CSRMATRIX_H
#define CSRMATRIX_H

#include "Matrix.h"

class CSRMatrix : public Matrix{
public:
    CSRMatrix(int rows, int columns, int *ind_ptr, int *indices, double *values_ptr);
    ~CSRMatrix();
};

#endif //CSRMATRIX_H
