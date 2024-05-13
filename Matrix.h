#pragma once
#include <vector>

class Matrix
{
public:
	int N;
	std::vector<std::vector<double>> matrix;
	Matrix operator * (double a);
	Matrix operator + (Matrix B);
	Matrix(int _N);
	Matrix();
};