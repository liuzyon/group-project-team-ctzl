#include <iostream>
#include<vector>
#include <cmath>
#include"Matrix.h"
#include"test.h"


Test::Test(/* args */)
{
}
Test::~Test() {

}


bool Test::test(Matrix& A, std::vector<double> b, std::vector<double> x) {
	//initialise output for right multiply
	std::vector<double> output;
	output.resize(b.size());
	
	bool ans = false;
	
	//do multiply A*x
	for (int i = 0; i < A.cols; i++) {
		for (int j = 0; j < A.cols; j++) {
			output[i] += A.values[i * A.cols + j] * x[j];
		}
	}
	//calculate errors and add tolerance
	auto* error = new double[b.size()];
	for (int i = 0; i < b.size(); i++) {
		error[i] = output[i] - b[i];
		if (abs(error[i]) < 1e-3) {
			ans = true;
		}
	}

	delete[] error;
	return ans;
}