#include <iostream>
#include "Matrix.h"
#include "DenseJacobi.h"
#include "DenseGaussSeidel.h"

int main()
{
    int rows = 5;
    int cols = 5;

    // Testing our matrix class
    auto *dense_mat = new Matrix(rows, cols, true);
    // Now we need to go and fill our matrix
    // Let's just fill it with random values for now
    for (int i = 0; i < rows * cols; i++)
    {
      dense_mat->values[i] = i;
    }
    // Now let's test printing our matrix
    dense_mat->printValues();
    std::cout << "A input:\n";
    // Now let's test printing our matrix with our other function
    dense_mat->printMatrix();


    std::vector<double> b = {1,2,3,4,5};
    std::cout << "b input:\n";
    for (int i = 0; i < b.size(); i++)
    {
        std::cout << b[i] << ' ';
    }
    DenseJacobi sv;
    std::vector<double> x = sv.solve(*dense_mat, b) ;

    
    std::cout << "\nx output:\n";
    for (int i = 0; i < x.size(); i++)
    {
       std::cout << x[i] << ' ';
    }

    // output for dense Gauss Seidel solver
    DenseGaussSeidel sv2;
    std::vector<double> x2 = sv2.solve(*dense_mat, b) ;
    
    std::cout << "\nx output:\n";
    for (int i = 0; i < x2.size(); i++)
    {
       std::cout << x2[i] << ' ';
    }


    // Don't forget to call delete (ie the destructor) once we're done, otherwise
    // we will leak memory from the new called inside our matrix
    delete dense_mat;
}