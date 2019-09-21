#include<iostream>
#include<math.h>
#include"matrix.h"

using namespace std;
class R {
public:
	R(float a, float b, float c, float d, float e, float f, float g, float h, float i) {
		a1 = a; a2 = b; a3 = c; b1 = d; b2 = e; b3 = f; c1 = g; c2 = h; c3 = i;
	}
	void print() {
		cout << "R矩阵为" << endl << a1 << " " << a2 << " " << a3 << endl << b1 << " " << b2 << " " << b3 << endl << c1 << " " << c2 << " " << c3 << endl;
	}
	float a1, a2, a3, b1, b2, b3, c1, c2, c3;
};

class PointPair {
public:
	PointPair(float x1, float y1, float X1, float Y1, float Z1) {
		x = x1/1000; y = y1/1000; X = X1; Y = Y1; Z = Z1;
	}
	Matrix getMatA(R R, Matrix mX, float f) {
		float Zba = R.a3 * (X - mX.data[0]) + R.b3 * (Y - mX.data[1]) + R.c3 * (Z - mX.data[2]);
		float omega = mX.data[4];
		float kappa = mX.data[5];
		float a11 = (R.a1 * f + R.a3 * x) / Zba;
		float a12 = (R.b1 * f + R.b3 * x) / Zba;
		float a13 = (R.c1 * f + R.c3 * x) / Zba;
		float a14 = y * sin(omega) - (x / f * (x * cos(kappa) - y * sin(kappa)) + f * cos(kappa)) * cos(omega);
		float a15 = -f * sin(kappa) - x / f * (x * sin(kappa) + y * cos(kappa));
		float a16 = y;
		float a21 = (R.a2 * f + R.a3 * y) / Zba;
		float a22 = (R.b2 * f + R.b3 * y) / Zba;
		float a23 = (R.c2 * f + R.c3 * y) / Zba;
		float a24 = -x * sin(omega) - (y / f * (x * cos(kappa) - y * sin(kappa)) - f * sin(kappa)) * cos(omega);
		float a25 = -f * cos(kappa) - y / f * (x * sin(kappa) + y * cos(kappa));
		float a26 = -x;
		float Result[] = { a11,a12,a13,a14,a15,a16,a21,a22,a23,a24,a25,a26 };
		Matrix r(Result, 2, 6);
		return r;
	}
	Matrix getMatL(R mT, Matrix mX, float f) {
		float XS, YS, ZS;
		XS = mX.data[0]; YS = mX.data[1]; ZS = mX.data[2];
		float x0 = -f * (mT.a1 * (X - XS) + mT.b1 * (Y - YS) + mT.c1 * (Z - ZS)) / (mT.a3 * (X - XS) + mT.b3 * (Y - YS) + mT.c3 * (Z - ZS));
		float y0 = -f * (mT.a2 * (X - XS) + mT.b2 * (Y - YS) + mT.c2 * (Z - ZS)) / (mT.a3 * (X - XS) + mT.b3 * (Y - YS) + mT.c3 * (Z - ZS));
		float l[] = { x - x0,y - y0 };
		Matrix L(l, 2, 1);
		return L;
	}
	float x, y;
	float X, Y, Z;
};

R getTranMat(float phi, float omega, float kappa) {
	float a1 = cos(phi) * cos(kappa) - sin(phi) * sin(omega) * sin(kappa);
	float a2 = -cos(phi) * sin(kappa) - sin(phi) * sin(omega) * cos(kappa);
	float a3 = -sin(phi) * cos(omega);
	float b1 = cos(omega) * sin(kappa);
	float b2 = cos(omega) * cos(kappa);
	float b3 = -sin(omega);
	float c1 = sin(phi) * cos(kappa) + cos(phi) * sin(omega) * sin(kappa);
	float c2 = -sin(phi) * sin(kappa) + cos(phi) * sin(omega) * cos(kappa);
	float c3 = cos(phi) * cos(omega);
	R a(a1, a2, a3, b1, b2, b3, c1, c2, c3);
	return a;
}

void solution(PointPair p1, PointPair p2, PointPair p3, PointPair p4, float f) {
	f /= 1000;
	//得到初始值
	float X = (p1.X + p2.X + p3.X + p4.X) / 4;
	float Y = (p1.Y + p2.Y + p3.Y + p4.Y) / 4;
	float Z = (p1.X - p2.X) / (p1.x - p2.x) * f;
	float phi = 0;
	float omega = 0;
	float kappa = 0;
	float x[] = { X, Y, Z, phi, omega, kappa };
	Matrix mX(x, 6, 1);
	cout << "初始值为" << endl;
	mX.transpose().printMat();
	for (int g = 0; g<10; g++){
		//计算像点矩阵L和每点的矩阵A
		R mT = getTranMat(mX.data[3], mX.data[4], mX.data[5]);
		Matrix L1 = p1.getMatL(mT, mX, f); Matrix A1 = p1.getMatA(mT, mX, f);
		Matrix L2 = p2.getMatL(mT, mX, f); Matrix A2 = p2.getMatA(mT, mX, f);
		Matrix L3 = p3.getMatL(mT, mX, f); Matrix A3 = p3.getMatA(mT, mX, f);
		Matrix L4 = p4.getMatL(mT, mX, f); Matrix A4 = p4.getMatA(mT, mX, f);
		Matrix L = appendRowMat(L1, L2, L3, L4);
		Matrix A = appendRowMat(A1, A2, A3, A4);
		//用delta=(ATA)-1(ATl)公式计算改正值
 		Matrix delta = multiplyMat(multiplyMat(A.transpose(),A).inverse(),multiplyMat(A.transpose(),L));
		mX.setMat(addMat(mX, delta));
		cout << "第" << g + 1 << "次迭代" << endl;
		mX.transpose().printMat();
		if (abs(delta.data[0])+ abs(delta.data[1]) + abs(delta.data[2]) < 0.5) {
			cout << "迭代结束" << endl << endl;
			getTranMat(mX.data[3], mX.data[4], mX.data[5]).print();
			return;
		}
	}
	cout << "迭代超过10次，程序已终止";
}

int main() {
	PointPair p1(-86.15, -68.99, 36589.41, 25273.32, 2195.17);
	PointPair p2(-53.4, 82.21, 37631.08, 31324.51, 728.69);
	PointPair p3(-14.78, -76.63, 39100.97, 24934.98, 2386.5);
	PointPair p4(10.46, 64.43, 40426.54, 30319.81, 757.31);
	solution(p1, p2, p3, p4, 153.24);
	getchar();
}
