#include "DenseGaussianElimination.h"
#include<iostream>
DenseGaussE::DenseGaussE(/* args */)
{
}

DenseGaussE::~DenseGaussE()
{
}

std::vector<double> DenseGaussE::solve(Matrix& A, std::vector<double>& b) {
	// This method is direct calculate the answer of Ax=b. 
	// This solver is only for the matix which exist the only answer x and no 0 on diagnal entries.
	// The computaional cost of Gaussian Elimination is O(n^3).
	//and we will also highlight issues that arise when the value is non-zero but small in pivoting.
	std::vector<double> x;
	x.resize(b.size());
	//make sure the sizes of matixs are what we can calculate
	if (A.cols != A.rows || A.cols != b.size()) {
		std::cerr << "Input dimensions for matrices don't match" << std::endl;
		return x;
	}
	//triangularise
	// covert A into upper triangluar form through row operations.
	//The same row operations are performed on the vector b.
	double s;
	// Loop over each pivot row all but the last row which we will never need to use as a pivot
	for (int k = 0; k < A.rows -1; k++) {
		//Loop over each row below the pivot row, including the last row which we do need to update
		for (int i = k + 1; i < A.rows; i++) {
			//Define the scaling factor for this row outside the innermost loop otherwise
			//There's also a performance saving from not recomputing things when not strictly necessary
			s = A.values[i * A.cols + k] / A.values[k * A.cols + k];
			//Update the current row of A by looping over the column j
			//start the loop from k as we can assume the entries before this are already zero
			for (int j = k; j < b.size(); j++) {
				A.values[i * A.cols + j] = A.values[i * A.cols + j] - s * A.values[k * A.cols + j];
			}
			// and update the corresponding entry of b
			b[i] = b[i] - s * b[k];
		}
	}

	//backsubstitude

	//start at the end (row n-1) and work backwards
	for (int k = b.size() - 1; k > -1; k--) {
		s = 0;
		for (int j = k + 1; j < b.size(); j++) {
			//sum of Akj*xj from k+1 to n
			s = s + A.values[k * A.cols + j] * x[j];
		}
		//xk = (1/akk)*(bk - sum of Akj*xj from k+1 to n)
		x[k] = (b[k] - s) / A.values[k * A.cols + k];
	}

	return x;
}