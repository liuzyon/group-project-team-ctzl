#ifndef TEST_H
#define TEST_H

#include"Matrix.h"
class Test 
{
private:
    /* data */
public:
    Test(/* args */);
    ~Test();
    bool test(Matrix& A, std::vector<double> b, std::vector<double> x);
};

#endif
