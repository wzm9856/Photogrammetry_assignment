#pragma once

class Matrix{
public:
	Matrix(float* mat, int row_size, int col_size);
	void setMat(Matrix mat);
	float* getMat();
	Matrix transpose();
	Matrix inverse();
	void printMat();
	int rowsize, colsize;
	float data[100];
};

void printMat(float* data, int rowsize, int colsize);
Matrix multiplyMat(Matrix mat1, Matrix mat2);
Matrix minusMat(Matrix mat1, Matrix mat2);
Matrix addMat(Matrix mat1, Matrix mat2);
Matrix appendRowMat(Matrix mat1, Matrix mat2, Matrix mat3, Matrix mat4);
