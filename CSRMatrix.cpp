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

//   // Check if our output matrix has had space allocated to it
//   if (output.values != nullptr)
//   {
//      // Check our dimensions match
//      if (this->rows != mat_left.cols || mat_left.rows != output.rows)
//      {
//         std::cerr << "Input dimensions for matrices don't match" << std::endl;
//         return;
//      }
//   }
//   // The output hasn't been preallocated, so we are going to do that
//   else
//   {
//      std::cerr << "OUTPUT HASN'T BEEN ALLOCATED" << std::endl;
//
//   }

   // HOW DO WE SET THE SPARSITY OF OUR OUTPUT MATRIX HERE??


    // Store values and its index in row and col
    std::vector<int> res_row_index;
    std::vector<int> res_col_index;
    std::vector<T> res_values;


    // Left：mat_left right：this
    // Iterate all values in mat_left
    for (int i = 0; i < mat_left.nnzs; ++i)
    {
        // For each value in left
        for (int j = 0; j < this->nnzs; ++j)
        {
            // Check if row index of the value in right matches the col index in left
            if (mat_left.col_index[i] == this->calRowIndex(j))
            {
                T value = mat_left.values[i] * this->values[j];
                int row = mat_left.calRowIndex(i);  // result row: the row index of left value
                int col = this->col_index[j];   // result column: the col index of right value

                // Iterate all values in result values list
                bool isExist = false;
                for (int k = 0; k < res_values.size(); ++k)
                {
                    // If it has value which has same row index and col index, merge the value.
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

    // delete 0 item in values and its index in col_index and row_index
    for (int i = 0; i < res_values.size(); ++i)
    {
        if (res_values[i] == 0) {
            res_values.erase(res_values.begin() + i);
            res_col_index.erase(res_col_index.begin() + i);
            res_row_index.erase(res_row_index.begin() + i);
        }
    }

    // construct values
    T *values = new T[res_values.size()];
    for (int i = 0; i < res_values.size(); ++i)
    {
        values[i] = res_values[i];
    }


    // construct col_index
    int *col_index = new int[res_col_index.size()];
    for (int i = 0; i < res_col_index.size(); ++i)
    {
        col_index[i] = res_col_index[i];
    }



    // construct row_position
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

}

// calculate the row index for the value belongs to
template <class T>
int CSRMatrix<T>::calRowIndex(int value_index)
{
    for (int i = 1; i < this->rows+1; ++i)
    {
        if (value_index < this->row_position[i]) {
            return i-1;
        }
    }
    return -1;
}
