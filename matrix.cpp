#include<iostream>
#include"matrix.h"

Matrix::Matrix(float* mat, int row_size, int col_size){
	if (row_size * col_size > 100)throw "Matrix is too large!";
	rowsize = row_size;
	colsize = col_size;
	memcpy(data, mat, rowsize * colsize * sizeof(float));
}

void Matrix::setMat(Matrix mat) {
	rowsize = mat.rowsize;
	colsize = mat.colsize;
	memcpy(data, mat.data, rowsize * colsize * sizeof(float));
}

float* Matrix::getMat(){
	return data;
}

Matrix Matrix::transpose() {
	float* Result = new float[colsize * rowsize];
	for (int i = 0; i < rowsize; i++) {
		for (int j = 0; j < colsize; j++) {
			Result[j * rowsize + i] = data[i * colsize + j];
		}
	}
	Matrix r(Result, colsize, rowsize);
	delete[] Result;
	return r;
}

//求逆矩阵
//我们使用LU分解求逆乘回去的方法
//来不及了就不允许行变换也没检查是不是可逆矩阵了。。。
Matrix Matrix::inverse() {
	if (rowsize != colsize) throw "Not a square matrix!";
	int size = rowsize;
	int size2 = size * size;

	float* L = new float[size2];	memset(L, 0, size2 * sizeof(float));
	float* U = new float[size2];	memcpy(U, data, size2 * sizeof(float));
	for (int i = 0; i < size; i++) {  //填充L第i列
		for (int j = i + 1; j < size; j++) {  //填充L(j,i),改变U第j行
			L[j * size + i] = U[j * size + i] / U[i * size + i];
			for (int k = i; k < size; k++) {   //改变U(j,k)
				U[j * size + k] -= U[i * size + k] * L[j * size + i];
			}
		}
	}
	for (int i = 0; i < size; i++) {
		L[i * size + i] = 1;
	}

	//对L矩阵求逆
	float* Li = new float[size2]; memset(Li, 0, size2 * sizeof(float));
	for (int i = 0; i < size; i++) {
		Li[i * size + i] = 1;
	}
	for (int i = 1; i < size; i++) {
		for (int j = i - 1; j >= 0; j--) {
			for (int k = j + 1; k <= i; k++) {
				Li[i * size + j] -= Li[i * size + k] * L[k * size + j];
			}
		}
	}

	//对U矩阵求逆
	float* Ui = new float[size2];
	memset(Ui, 0, size2 * sizeof(float));
	for (int i = 0; i < size; i++) {
		Ui[i * size + i] = 1 / U[i * size + i];
	}
	for (int i = size - 2; i >= 0; i--) {
		for (int j = i + 1; j < size; j++) {
			for (int k = j - 1; k >= i; k--) {
				Ui[i * size + j] -= Ui[i * size + k] * U[k * size + j];
			}
			Ui[i * size + j] /= U[j * size + j];
		}
	}

	//Ui*Li
	float* Result = new float[size2];
	memset(Result, 0, size2 * sizeof(float));
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			for (int k = 0; k < size; k++) {
				Result[i * size + j] += Ui[i * size + k] * Li[k * size + j];
			}
		}
	}

	delete[] L, U, Li, Ui;
	Matrix r(Result, rowsize, colsize);
	return r;
}

//测试用函数 方便打印矩阵
void Matrix::printMat() {
	for (int i = 0; i < rowsize; i++) {
		for (int j = 0; j < colsize; j++) {
			std::cout << data[i * colsize + j] << '\t';
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
}

//直接打印内部数据
void printMat(float* data, int rowsize, int colsize) {
	for (int i = 0; i < rowsize; i++) {
		for (int j = 0; j < colsize; j++) {
			std::cout << data[i * colsize + j] << " ";
		}
		std::cout << std::endl;
	}
}

//矩阵相乘
Matrix multiplyMat(Matrix mat1, Matrix mat2) {
	if (mat1.colsize != mat2.rowsize) throw "Can't multiply these two matrix!";	 
	int rowsize = mat1.rowsize;
	int colsize = mat2.colsize;
	int emm = mat1.colsize;
	float* Result = new float[rowsize * colsize];
	memset(Result, 0, rowsize * colsize * sizeof(float));
	for (int i = 0; i < rowsize; i++) {
		for (int j = 0; j < colsize; j++) {
			for (int k = 0; k < emm; k++) {
				Result[i * colsize + j] += mat1.data[i * mat1.colsize + k] * mat2.data[k * mat2.colsize + j];
			}
		}
	}
	Matrix r(Result, rowsize, colsize);
	delete[] Result;
	return r;
}

//矩阵相减
Matrix minusMat(Matrix mat1, Matrix mat2) {
	if (mat1.colsize != mat2.colsize || mat1.rowsize != mat2.rowsize)
		throw "Can't minus these two matrix";
	int rowsize = mat1.rowsize;
	int colsize = mat2.colsize;
	float* Result = new float[rowsize * colsize];
	for (int i = 0; i < rowsize; i++) {
		for (int j = 0; j < colsize; j++) {
			Result[i * colsize + j] = mat1.data[i * colsize + j] - mat2.data[i * colsize + j];
		}
	}
	Matrix r(Result, rowsize, colsize);
	delete[] Result;
	return r;
}

//矩阵相加
Matrix addMat(Matrix mat1, Matrix mat2) {
	if (mat1.colsize != mat2.colsize || mat1.rowsize != mat2.rowsize)
		throw "Can't add these two matrix";
	int rowsize = mat1.rowsize;
	int colsize = mat2.colsize;
	float* Result = new float[rowsize * colsize];
	for (int i = 0; i < rowsize; i++) {
		for (int j = 0; j < colsize; j++) {
			Result[i * colsize + j] = mat1.data[i * colsize + j] + mat2.data[i * colsize + j];
		}
	}
	Matrix r(Result, rowsize, colsize);
	delete[] Result;
	return r;
}

//叠加矩阵
Matrix appendRowMat(Matrix mat1, Matrix mat2, Matrix mat3, Matrix mat4) {
	if (mat1.colsize != mat2.colsize || mat2.colsize != mat3.colsize || mat3.colsize != mat4.colsize) 
		throw "Can't append these matrixes!";
	int rowsize = mat1.rowsize + mat2.rowsize + mat3.rowsize + mat4.rowsize;
	int colsize = mat1.colsize;
	float* Result = new float[rowsize * colsize];
	memcpy(Result, mat1.data, mat1.rowsize * colsize * sizeof(float));
	memcpy(Result + mat1.rowsize * colsize, mat2.data, mat2.rowsize * colsize * sizeof(float));
	memcpy(Result + (mat1.rowsize + mat2.rowsize) * colsize, mat3.data, mat3.rowsize * colsize * sizeof(float));
	memcpy(Result + (mat1.rowsize + mat2.rowsize + mat3.rowsize) * colsize, mat4.data, mat4.rowsize * colsize * sizeof(float));
	Matrix r(Result, rowsize, colsize);
	delete[] Result;
	return r;
}