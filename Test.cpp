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
	// use the A*x_-b = r
	//caluculate the norm of r, if the norm of r is less then 1e-5
	// this x_ pass our test and return 1.


	//initialise output for right multiply
	std::vector<double> output;
	output.resize(b.size());
	
	bool ans = false;
	
	//do multiply A*x_
	for (int i = 0; i < A.cols; i++) {
		for (int j = 0; j < A.cols; j++) {
			output[i] += A.values[i * A.cols + j] * x[j];
		}
	}
	//calculate error 
	double s = 0;
	for (int i = 0; i < b.size(); i++) {
		s += pow((output[i] - b[i]),2);
	}
	if (sqrt(s) < 1e-5) {
		ans = true;
	}

	
	return ans;
}