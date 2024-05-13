#include "Matrix.h"

Matrix::Matrix(int _N)
{
	matrix.resize(_N);
	for (int i = 0; i < _N; i++)
		matrix[i].resize(_N);
	N = _N;
}
Matrix::Matrix()
{
	matrix.resize(4);
	for (int i = 0; i < 4; i++)
		matrix[i].resize(4);
	N = 4;
}
Matrix Matrix::operator+(Matrix B)
{
	if (B.N != N)
		return -1;
	Matrix C(N);

	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
			C.matrix[i][j] = matrix[i][j] + B.matrix[i][j];

	return C;
}
Matrix Matrix::operator*(double a)
{
	Matrix C(N);
	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
			C.matrix[i][j] = a * matrix[i][j];
	return C;
}
