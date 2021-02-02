#include <iostream>
#include <vector>
#include "Matrix.h"
#include "Matrix.cpp"
#include "CSRMatrix.h"
#include "CSRMatrix.cpp"
#include "Solver.h"
#include "Solver.cpp"

using namespace std;

int main()
{
    // please do not open the first and second tests at the same time.
    //welcome to  use test method to your solver.
    //~~~~~~~~~~~~~~~~~
    //~~~~test one~~~~~
    //~~~~~~~~~~~~~~~~~
    //specific test for DenseGaussianEliminationpp solver
//    int rows = 3;
//    int cols = 3;
//    double input[9] = { 2, 3, -4, 6, 8, 2, 4, 8, -6 };
//    auto* A_GEpp = new Matrix(rows, cols, input);
//    std::cout << "A input:\n";
//    A_GEpp->printMatrix();
//    std::vector<double> b = { 5., 3., 19. };
//    std::cout << "b input:\n";
//    for (int i = 0; i < b.size(); i++)
//    {
//        std::cout << b[i] << ' ';
//    }
//
//    DenseGaussEPP gepp;
//    std::vector<double> xgepp = gepp.solve(*A_GEpp, b);
//
//    std::cout << "\nx output:\n";
//    for (int i = 0; i < xgepp.size(); i++)
//    {
//        std::cout << xgepp[i] << ' ';
//    }
//
//    Test t;
//    bool res = t.test(*A_GEpp, b, xgepp);
//    std::cout << "\n" << res;
//    delete A_GEpp;
    //~~~~~~~~~~~~~~~~~
    //~~~~test two~~~~~
    //~~~~~~~~~~~~~~~~~
    //specific test for DenseGaussianEliminationpp solver
    /*int rows = 3;
    int cols = 3;
    double input[9] = { 2, 3, -4, 3, -1, 2, 4, 2, 2 };
    auto* A_GE = new Matrix(rows, cols, input);
    std::cout << "A input:\n";
    A_GE->printMatrix();
    std::vector<double> b = { 10., 3., 8. };
    std::cout << "b input:\n";
    for (int i = 0; i < b.size(); i++)
    {
        std::cout << b[i] << ' ';
    }

    DenseGaussEPP ge;
    std::vector<double> xge = ge.solve(*A_GE, b);
    std::cout << "\nx output:\n";
    for (int i = 0; i < xge.size(); i++)
    {
        std::cout << xge[i] << ' ';
    }
    Test t;
    bool res = t.test(*A_GE, b, xge);
    std::cout << "\n" << res;
    delete A_GE;*/
    //~~~~~~~~~~~~~~~~~
    //~~~~test end~~~~~
    //~~~~~~~~~~~~~~~~~
    //int rows = 5;
    //int cols = 5;

    //// Testing our matrix class
    //auto *dense_mat = new Matrix(rows, cols, true);
    //// Now we need to go and fill our matrix
    //// Let's just fill it with random values for now
    //for (int i = 0; i < rows * cols; i++)
    //{
    //  dense_mat->values[i] = i;
    //}
    //// Now let's test printing our matrix
    //dense_mat->printValues();
    //std::cout << "A input:\n";
    //// Now let's test printing our matrix with our other function
    //dense_mat->printMatrix();


    //std::vector<double> b = {1,2,3,4,5};
    //std::cout << "b input:\n";
    //for (int i = 0; i < b.size(); i++)
    //{
    //    std::cout << b[i] << ' ';
    //}
    //DenseJacobi sv;
    //std::vector<double> x = sv.solve(*dense_mat, b) ;

    //
    //std::cout << "\nx output:\n";
    //for (int i = 0; i < x.size(); i++)
    //{
    //   std::cout << x[i] << ' ';
    //}

    //// output for dense Gauss Seidel solver
    //DenseGaussSeidel sv2;
    //std::vector<double> x2 = sv2.solve(*dense_mat, b) ;
    //
    //std::cout << "\nx output:\n";
    //for (int i = 0; i < x2.size(); i++)
    //{
    //   std::cout << x2[i] << ' ';
    //}

    //delete dense_mat;

    ////CSRMatrix
    //int rows = 3;
    //int columns = 3;
    //int a[4] = {0, 2, 3, 6};
    //int b[6] = {0, 2, 2, 0, 1, 2};
    //double c[6] = {1, 2, 3, 4, 5, 6};
    //int *int_ptr = a;
    //int *indices = b;
    //double *values_ptr= c;
    //auto *CSR_mat = new CSRMatrix(rows, columns, int_ptr, indices, values_ptr);
    //CSR_mat->printMatrix();
    //delete CSR_mat;

    // Test LUFactorisation
//    auto* A_LUFactorisation = new Matrix(rows, cols, input);
//    std::vector<double> b = { 2., 1., 5. };
//    Solver<double> sv;
//    std::vector<double> x = sv.DenseLUFactorisationSolve(*A_LUFactorisation, b);
//    // expect output: -4.375 -2.75 2.5
//    for (int i = 0; i < x.size(); i++)
//    {
//        std::cout << x[i] << ' ';
//    }

    cout << "Dense Matrix" << endl;
    int rows = 3;
    int cols = 3;
    double input[9] = { 4, -3, 1, -2, 3, -2, 4, 0, 9};
    Matrix<double> A_dense(rows, cols, input);
    double *b = new double [3]{ 2., 1., 5.};
    double *x = new double [3];
    Solver<double> sv;
    sv.DenseLUFactorisationSolve(A_dense, b, x);
    for (int i = 0; i < 3; ++i)
    {
        cout << x[i] << " ";
    }

    cout << endl;
    // if no input for user tolerance, default is 1e-6
    sv.dense_jacobi_solver(A_dense, b, x, 0.01);
    for (int i = 0; i < 3; ++i)
    {
        cout << x[i] << " ";
    }

    cout << endl;
    // if no input for user tolerance, default is 1e-6
    sv.dense_gauss_seidel_solver(A_dense, b, x, 0.1);
    for (int i = 0; i < 3; ++i)
    {
        cout << x[i] << " ";
    }

    delete[] b;
    delete[] x;

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // test for sparse mmatrix
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    cout << endl << endl;
    int sparse_rows = 5;
    int sparse_cols = 5;
    int nnzs = 12;

    double sparse_values[12] = {10,2,3,40,5,6,70,8,9,10,110,120};
    int row_position[6] = {0,2,5,9,11,12};
    int col_index[12] = {0,3,0,1,3,0,2,3,4,2,3,4};
   
    auto *sparse_mat = new CSRMatrix<double>(sparse_rows, sparse_cols, nnzs, sparse_values, row_position, col_index);

    // Now let's print it
    sparse_mat->printMatrix();
    cout << endl;


    double* sparse_b = new double[5]{5,3,8,1,4};
    double* sparse_x = new double[5];
    // store sparse matvec result
    double* sparse_output = new double[5];


    // if no input for user tolerance, default is 1e-6
    sv.sparse_jacobi_solver(*sparse_mat, sparse_b, sparse_x, 0.0001);

    cout << "Sparse Jacobi result (tol = 1e-4): " << endl;
    for (int i = 0; i < sparse_cols; i++)
    {
       cout << sparse_x[i] << " ";
    }
    cout << endl << endl;
    sparse_mat->matVecMult(sparse_x, sparse_output);
    // Let's print the results
    cout << "Mat vec using above result " << endl;
    for (int i = 0; i < sparse_rows; i++)
    {
       cout << " " << sparse_output[i];
    }
    cout << endl << endl;


   // if no input for user tolerance, default is 1e-6
   sv.sparse_gauss_seidel_solver(*sparse_mat, sparse_b, sparse_x);
   cout << "Sparse Gauss Seidel result: (tol = 1e-6)" << endl;
   for (int i = 0; i < sparse_cols; i++)
   {
      cout << sparse_x[i] << " ";
   }
   cout << endl << endl;
   sparse_mat->matVecMult(sparse_x, sparse_output);
   // Let's print the results
   cout << "Mat vec using above result " << endl;
   for (int i = 0; i < sparse_rows; i++)
   {
      cout << " " << sparse_output[i];
   }
   cout << endl;

   delete sparse_mat;
   delete[] sparse_output;
   delete[] sparse_x;
   delete[] sparse_b;

}