#include "DenseGaussianEliminationpp.h"
#include<iostream>
#include<cmath>
DenseGaussEPP::DenseGaussEPP(/* args */)
{
}

DenseGaussEPP::~DenseGaussEPP()
{
}

std::vector<double> DenseGaussEPP::solve(Matrix& A, std::vector<double>& b) {
	// This method is direct calculate the answer of Ax=b. 
	// This solver is only for the matix which exist the only answer x.
	// The computaional cost of Gaussian Elimination is O(n^3).
	//and we fix the issue we meet in DenseGaussE class we use partial pivoting.
	//This is an upgrade of Gaussian Elimination method.
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
	int kmax;
	// Loop over each pivot row all but the last row which we will never need to use as a pivot
	for (int k = 0; k < A.rows - 1; k++) {
		//Swap rows so we are always dividing through by the largest magnitude number.
		//initiatise kmax with the current pivot row (k)
		kmax = k;
		//loop over all entries below the pivot and select the k with the largest abs value
		for (int i = k + 1; i < A.rows; i++) {
			if (abs(A.values[kmax * A.cols + k]) < abs(A.values[i * A.cols + k])) {
				kmax = i;
			}
				
		}
		//and swap the current pivot row (k) with the row with the largest abs value below the pivot
		std::vector<double> iA;
		iA.resize(A.cols);
		for (int j = 0; j < A.cols; j++) {
			iA[j] = A.values[kmax * A.cols + j];
		}
		double ib = b[kmax];

		for (int j = 0; j < A.cols; j++) {
			A.values[kmax * A.cols + j] = A.values[k * A.cols + j];
		}
		b[kmax] = b[k];
		for (int j = 0; j < A.cols; j++) {
			A.values[k * A.cols + j] = iA[j];
		}
		b[k] = ib;


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