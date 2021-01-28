//
// Created by Zhiyong Liu on 2021/1/27.
//

/*#include "CSRMatrix.h"

CSRMatrix::CSRMatrix(int rows, int columns, int *ind_ptr, int *indices, double *values_ptr)
{

    this->rows = rows;
    this->cols = columns;
    this->size_of_values = rows * columns;
    this->values = new double[size_of_values]();

    for (int i = 0; i < rows; ++i)
    {
        // i 为行index
        int row_start = ind_ptr[i];
        int row_end = ind_ptr[i+1];
        for (int j = row_start; j < row_end; ++j)
        {
            // indices[j] 为对应元素所对应行偏移量
            values[i*columns + indices[j]] = values_ptr[j];
        }
    }
}

CSRMatrix::~CSRMatrix()
{
    delete[] this->values;
}*/
