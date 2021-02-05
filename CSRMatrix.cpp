#include <iostream>
#include<vector>
#include "CSRMatrix.h"

// Constructor - using an initialisation list here
template <class T>
CSRMatrix<T>::CSRMatrix(int rows, int cols, int nnzs, bool preallocate): Matrix<T>(rows, cols, false), nnzs(nnzs)
{
   // If we don't pass false in the initialisation list base constructor, it would allocate values to be of size
   // rows * cols in our base matrix class
   // So then we need to set it to the real value we had passed in
   this->preallocated = preallocate;

   // If we want to handle memory ourselves
   if (this->preallocated)
   {
      // Must remember to delete this in the destructor
      this->values = new T[this->nnzs];
      this->row_position = new int[this->rows+1];
      this->col_index = new int[this->nnzs];
   }
}

// Constructor - now just setting the value of our T pointer
template <class T>
CSRMatrix<T>::CSRMatrix(int rows, int cols, int nnzs, T *values_ptr, int *row_position, int *col_index): Matrix<T>(rows, cols, values_ptr), nnzs(nnzs), row_position(row_position), col_index(col_index)
{}

// destructor
template <class T>
CSRMatrix<T>::~CSRMatrix()
{
   // Delete the values array
   if (this->preallocated){
      delete[] this->row_position;
      delete[] this->col_index;
   }
   // The super destructor is called after we finish here
   // This will delete this->values if preallocated is true
}

// Explicitly print out the values in values array as if they are a matrix
template <class T>
void CSRMatrix<T>::printMatrix() 
{ 
   std::cout << "Printing matrix" << std::endl;
   std::cout << "Values: ";
   for (int j = 0; j< this->nnzs; j++)
   {  
      std::cout << this->values[j] << " ";      
   }
   std::cout << std::endl;
   std::cout << "row_position: ";
   for (int j = 0; j< this->rows+1; j++)
   {  
      std::cout << this->row_position[j] << " ";      
   }
   std::cout << std::endl;   
   std::cout << "col_index: ";
   for (int j = 0; j< this->nnzs; j++)
   {  
      std::cout << this->col_index[j] << " ";      
   }
   std::cout << std::endl;   
}

// Do a matrix-vector product
// output = this * input
template<class T>
void CSRMatrix<T>::matVecMult(T *input, T *output)
{
   if (input == nullptr || output == nullptr)
   {
      std::cerr << "Input or output haven't been created" << std::endl;
      return;
   }

   // Set the output to zero
   for (int i = 0; i < this->rows; i++)
   {
      output[i] = 0.0;
   }

   // Loop over each row
   for (int i = 0; i < this->rows; i++)
   {
      // Loop over all the entries in this col
      for (int val_index = this->row_position[i]; val_index < this->row_position[i+1]; val_index++)
      {
         // This is an example of indirect addressing
         // Can make it harder for the compiler to vectorise!
         output[i] += this->values[val_index] * input[this->col_index[val_index]];

      }
   }
}


// Do matrix matrix multiplication
// output = mat_left * this
template <class T>
void CSRMatrix<T>::matMatMult(CSRMatrix<T>& mat_left, CSRMatrix<T>& output)
{

   // Check our dimensions match
   if (this->cols != output.cols)
   {
      std::cerr << "Input dimensions for matrices don't match" << std::endl;
      return;
   }

   // Check if our output matrix has had space allocated to it
   if (output.values != nullptr) 
   {
      // Check our dimensions match
      if (this->rows != mat_left.cols || mat_left.rows != output.rows)
      {
         std::cerr << "Input dimensions for matrices don't match" << std::endl;
         return;
      }      
   }
   // The output hasn't been preallocated, so we are going to do that
   else
   {
      std::cerr << "OUTPUT HASN'T BEEN ALLOCATED" << std::endl;

   }



   // HOW DO WE SET THE SPARSITY OF OUR OUTPUT MATRIX HERE??

    // 存储结果的三元组集
    std::vector<T> res_row_index;
    std::vector<T> res_col_index;
    std::vector<T> res_values;

//    std::tuple<int, int, T> res;



    // 左矩阵：mat_left 右矩阵：this
    // 遍历mat_left所有的元素(col_index和values的index是一一对应的)
    for (int i = 0; i < mat_left.nnzs; ++i)
    {
        // 对于左矩阵每个元素，用其列值与所有右矩阵元素的行值匹配
        for (int j = 0; j < this->nnzs; ++j)
        {
            // 如果第二个矩阵有元素可以和第一个矩阵该元素满足想乘条件
            if (mat_left.col_index[i] == this->calRowIndex(j))
            {
                T value = mat_left.values[i] * this->values[j];
                int row = mat_left.calRowIndex(i);  // 行为第一个矩阵元素的行值
                int col = this->col_index[j];   // 列为第二个矩阵元素的列值

                // 对已有的三元组进行遍历，查看新生成的是否可以合并

                bool isExist = false;
                for (int k = 0; k < res_values.size(); ++k)
                {
                    if (res_row_index[k] == row && res_col_index[k] == col) {
                        res_values[k] += value;
                        isExist = true;
                        break;
                    }
                }

                if (!isExist) {
                    res_row_index.push_back(row);
                    res_col_index.push_back(col);
                    res_values.push_back(value);
                }
            }
        }
    }

    // 删除结果中为0的元素
    for (int i = 0; i < res_values.size(); ++i)
    {
        if (res_values[i] == 0) {
            res_values.erase(res_values.begin() + i);
            res_col_index.erase(res_col_index.begin() + i);
            res_row_index.erase(res_row_index.begin() + i);
        }
    }

    // 构成values
    T *values = new T[res_values.size()];
    for (int i = 0; i < res_values.size(); ++i)
    {
        values[i] = res_values[i];
    }


    // 构成col_index
    int *col_index = new int[res_col_index.size()];
    for (int i = 0; i < res_col_index.size(); ++i)
    {
        col_index[i] = res_col_index[i];
    }



    // 构成row_position
    int *row_position = new int[mat_left.rows+1];
    row_position[0] = 0;
    for (int i = 0; i < mat_left.rows; ++i)
    {
        int times = std::count(res_row_index.begin(), res_row_index.end(), i);
        row_position[i+1] = row_position[i] + times;
    }




    output.nnzs = res_values.size();
    output.rows = mat_left.rows;
    output.cols = this->cols;

    output.values = values;
    output.col_index = col_index;
    output.row_position = row_position;

    // 打印结果value，用于debug，可注释掉
    std::cout << "values: " << std::endl;
    for (int i = 0; i < res_values.size(); ++i)
    {
        std::cout << values[i] << " ";
    }
    std::cout << std::endl;

    // 打印结果col_index，用于debug，可注释掉
    std::cout << "col values: " << std::endl;
    for (int i = 0; i < res_col_index.size(); ++i)
    {
        std::cout << col_index[i] << " ";
    }
    std::cout << std::endl;

    // 打印结果row position，用于debug，可注释掉
    std::cout << "row position: " << std::endl;
    for (int i = 0; i < mat_left.rows+1; ++i)
    {
        std::cout << row_position[i] << " ";
    }
    std::cout << std::endl;

}

// row_size: 行数，不是row_position的size
template <class T>
int CSRMatrix<T>::calRowIndex(int value_index)
{

    for (int i = 1; i < this->rows+1; ++i)
    {
        if (value_index < this->row_position[i]) {
            return i-1; // 返回所属行数(0:第一行，1：第二行。。。。)
        }
    }
}
