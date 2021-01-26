cmake_minimum_required(VERSION 3.17)
project(assessment1)

set(CMAKE_CXX_STANDARD 17)

add_executable(assessment1 main.cpp Jacobi.h Jacobi.cpp Matrix.h Matrix.cpp Solver.h Solver.cpp DenseJacobi.h DenseJacobi.cpp)