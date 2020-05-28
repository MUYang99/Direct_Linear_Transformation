// DLT.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

#include "pch.h"
#include <iostream>
#include<stdio.h>
#include"string.h"
#include<stdlib.h>
#include<fstream>
#include<vector>
#include<algorithm>
#include<iomanip>
#include <sstream> 
#include <string>
#include<math.h>

using namespace std;
#define pixsize 0.0051966
#pragma warning(disable:4996);

struct PointData {
	int PointsNumber;
	double x;
	double y;
	double X;
	double Y;
	double Z;
	double xl;
	double yl;
	double xr;
	double yr;
};

/****************************************矩阵运算函数*************************************/
//矩阵转置函数
void MatrixTranposition(double *MatrixOrigin, double *MatrixFinal, int m, int n) {
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			MatrixFinal[i*m + j] = MatrixOrigin[j*n + i];
		}
	}
}

//矩阵相乘
void MatrixMultiplication(double *MatrixOrign, double *MatrixOrignCopy, double *MatrixPlus, int m, int n, int s) {
	int i, j, k;
	for (i = 0; i < m; i++)
		for (j = 0; j < s; j++)
		{
			MatrixPlus[i*s + j] = 0.0;
			for (k = 0; k < n; k++)
				MatrixPlus[i*s + j] += MatrixOrign[i*n + k] * MatrixOrignCopy[j + k * s];
		}
}

//矩阵求逆函数
void MatrixInversion(double *MatrixInverse, int n) {
	int i, j, k;
	for (k = 0; k < n; k++)
	{
		for (i = 0; i < n; i++)
		{
			if (i != k)
				MatrixInverse[i*n + k] = -MatrixInverse[i*n + k] / MatrixInverse[k*n + k];
		}
		MatrixInverse[k*n + k] = 1 / MatrixInverse[k*n + k];
		for (i = 0; i < n; i++)
		{
			if (i != k)
			{
				for (j = 0; j < n; j++)
				{
					if (j != k)
						MatrixInverse[i*n + j] += MatrixInverse[k*n + j] * MatrixInverse[i*n + k];
				}
			}
		}
		for (j = 0; j < n; j++)
		{
			if (j != k)
				MatrixInverse[k*n + j] *= MatrixInverse[k*n + k];
		}
	}
}

/*************************************DLT运算函数*****************************************/
//Li的近似值的求解
void LiApproximation(double *&L, PointData *Locate, int n) { //Locate为物方控制点
	double *M = new double[n * 2 * 11];
	double *W = new double[n * 2];
	double *MT = new double[11 * 2 * n];
	double *MTM = new double[11 * 11];//由误差方程式得到
	double Result[11] = { 0.0 };
	for (int i = 0; i < n; i++) {
		M[i * 2 * 11 + 0] = Locate[i].X;
		M[i * 2 * 11 + 1] = Locate[i].Y;
		M[i * 2 * 11 + 2] = Locate[i].Z;
		M[i * 2 * 11 + 3] = 1.0;
		M[i * 2 * 11 + 4] = 0.0;
		M[i * 2 * 11 + 5] = 0.0;
		M[i * 2 * 11 + 6] = 0.0;
		M[i * 2 * 11 + 7] = 0.0;
		M[i * 2 * 11 + 8] = Locate[i].x*Locate[i].X;
		M[i * 2 * 11 + 9] = Locate[i].x*Locate[i].Y;
		M[i * 2 * 11 + 10] = Locate[i].x*Locate[i].Z;
		M[i * 2 * 11 + 11] = 0.0;
		M[i * 2 * 11 + 12] = 0.0;
		M[i * 2 * 11 + 13] = 0.0;
		M[i * 2 * 11 + 14] = 0.0;
		M[i * 2 * 11 + 15] = Locate[i].X;
		M[i * 2 * 11 + 16] = Locate[i].Y;
		M[i * 2 * 11 + 17] = Locate[i].Z;
		M[i * 2 * 11 + 18] = 1.0;
		M[i * 2 * 11 + 19] = Locate[i].y*Locate[i].X;
		M[i * 2 * 11 + 20] = Locate[i].y*Locate[i].Y;
		M[i * 2 * 11 + 21] = Locate[i].y*Locate[i].Z;//M矩阵的存储


		W[i * 2] = Locate[i].x;
		W[i * 2 + 1] = Locate[i].y;//W矩阵的存储

		MatrixTranposition(M, MT, 2 * n, 11);
		MatrixMultiplication(MT, M, MTM, 11, 2 * n, 11);
		MatrixInversion(MTM, 11);
		MatrixMultiplication(MT, W, Result, 11, 2 * n, 1);
		MatrixMultiplication(MTM, Result, L, 11, 11, 1);
	}
}

//内方位元素近似值求解
void InteriorElement(double &x0, double &y0, double &fx, double *L) {
	x0 = 0.0;
	y0 = 0.0;
	x0 = -(L[1] * L[9] + L[2] * L[10] + L[3] * L[11]) / (L[9] * L[9] + L[10] * L[10] + L[11] * L[11]);
	y0 = -(L[5] * L[9] + L[6] * L[10] + L[7] * L[11]) / (L[9] * L[9] + L[10] * L[10] + L[11] * L[11]);
	double m = 0.0;
	double n = 0.0;
	double s = 0.0;
	m = (L[0] * L[0] + L[1] * L[1] + L[2] * L[2]) / (L[8] * L[8] + L[9] * L[9] + L[10] * L[10]) - x0 * x0;
	n = (L[4] * L[4] + L[5] * L[5] + L[6] * L[6]) / (L[8] * L[8] + L[9] * L[9] + L[10] * L[10]) - y0 * y0;
	s = (L[0] * L[4] + L[1] * L[5] + L[2] * L[6]) / (L[8] * L[8] + L[9] * L[9] + L[10] * L[10]) - x0 * y0;
	fx = sqrt((m*n - s * s) / n);
}

//Li的精确值的求解以及内方位元素精确值的求解
void LiAccuration(double *&vv, double *&MTM, double *&LX, double &x0Accurate, double &y0Accurate, double &fAccurate, double *L, PointData *Locate, double x0, double y0, double f, int n) {
	double *A = new double[n];
	double *M = new double[2 * n * 13];
	double *W = new double[2 * n];
	vv = new double[2 * n];
	double *MT = new double[13 * 2 * n];
	MTM = new double[13 * 13];
	double *MTW = new double[13 * 1];
	double *ML = new double[2 * n];
	double r = 0.0;

	LX[0] = L[0];
	LX[1] = L[1];
	LX[2] = L[2];
	LX[3] = L[3];
	LX[4] = L[4];
	LX[5] = L[5];
	LX[6] = L[6];
	LX[7] = L[7];
	LX[8] = L[8];
	LX[9] = L[9];
	LX[10] = L[10];
	LX[11] = 0.0;//k1
	LX[12] = 0.0;//k2

	double mm = 0.0;
	double nn = 0.0;
	double ss = 0.0;
	double f_t = 0.0;
	double x0_t = 0.0;
	double y0_t = 0.0;
	x0Accurate = x0;
	y0Accurate = y0;
	double f_0 = 0.0;
	do {
		x0_t = x0Accurate;
		y0_t = y0Accurate;
		f_t = f_0;
		for (int i = 0; i < n; i++)
		{
			r = sqrt((Locate[i].x - x0Accurate)*(Locate[i].x - x0Accurate) + (Locate[i].y - y0Accurate)*(Locate[i].y - y0Accurate));//计算像点像径
			A[i] = LX[8] * Locate[i].X + LX[9] * Locate[i].Y + LX[10] * Locate[i].Z + 1;//计算A
			//计算M矩阵																
			M[i * 26] = -Locate[i].X / A[i];
			M[i * 26 + 1] = -Locate[i].Y / A[i];
			M[i * 26 + 2] = -Locate[i].Z / A[i];
			M[i * 26 + 3] = -1 / A[i];
			M[i * 26 + 4] = 0.0;
			M[i * 26 + 5] = 0.0;
			M[i * 26 + 6] = 0.0;
			M[i * 26 + 7] = 0.0;
			M[i * 26 + 8] = -Locate[i].x*Locate[i].X / A[i];
			M[i * 26 + 9] = -Locate[i].x*Locate[i].Y / A[i];
			M[i * 26 + 10] = -Locate[i].x*Locate[i].Z / A[i];
			M[i * 26 + 11] = -(Locate[i].x - x0Accurate)*r*r;
			M[i * 26 + 12] = -(Locate[i].x - x0Accurate)*r*r*r*r;
			M[i * 26 + 13] = 0.0;
			M[i * 26 + 14] = 0.0;
			M[i * 26 + 15] = 0.0;
			M[i * 26 + 16] = 0.0;
			M[i * 26 + 17] = -Locate[i].X / A[i];
			M[i * 26 + 18] = -Locate[i].Y / A[i];
			M[i * 26 + 19] = -Locate[i].Z / A[i];
			M[i * 26 + 20] = -1 / A[i];
			M[i * 26 + 21] = -Locate[i].y*Locate[i].X / A[i];
			M[i * 26 + 22] = -Locate[i].y*Locate[i].Y / A[i];
			M[i * 26 + 23] = -Locate[i].y*Locate[i].Z / A[i];
			M[i * 26 + 24] = -(Locate[i].y - y0Accurate)*r*r;
			M[i * 26 + 25] = -(Locate[i].y - y0Accurate)*r*r*r*r;
			//计算W矩阵
			W[i * 2] = Locate[i].x / A[i];
			W[i * 2 + 1] = Locate[i].y / A[i];
		}
		//计算法方程式
		MatrixTranposition(M, MT, 2 * n, 13);
		MatrixMultiplication(MT, M, MTM, 13, 2 * n, 13);
		MatrixInversion(MTM, 13);
		MatrixMultiplication(MT, W, MTW, 13, 2 * n, 1);
		MatrixMultiplication(MTM, MTW, LX, 13, 13, 1);
		//利用精确的L系数，计算内方位元素的精确值
		x0Accurate = -(LX[0] * LX[8] + LX[1] * LX[9] + LX[2] * LX[10]) / (LX[8] * LX[8] + LX[9] * LX[9] + LX[10] * LX[10]);
		y0Accurate = -(LX[4] * LX[8] + LX[5] * LX[9] + LX[6] * LX[10]) / (LX[8] * LX[8] + LX[9] * LX[9] + LX[10] * LX[10]);
		mm = (LX[0] * LX[0] + LX[1] * LX[1] + LX[2] * LX[2]) / (LX[8] * LX[8] + LX[9] * LX[9] + LX[10] * LX[10]) - x0Accurate * x0Accurate;
		nn = (LX[4] * LX[4] + LX[5] * LX[5] + LX[6] * LX[6]) / (LX[8] * LX[8] + LX[9] * LX[9] + LX[10] * LX[10]) - y0Accurate * y0Accurate;
		ss = (LX[0] * LX[4] + LX[1] * LX[5] + LX[2] * LX[6]) / (LX[8] * LX[8] + LX[9] * LX[9] + LX[10] * LX[10]) - x0Accurate * y0Accurate;
		fAccurate = sqrt((mm*nn - ss * ss) / nn);
		f_0 = fAccurate;
	} while (fabs(fAccurate - f_t) > 0.01 || fabs(x0Accurate - x0_t) > 0.01 || fabs(y0Accurate - y0_t) > 0.01);
	//计算误差方程式
	MatrixMultiplication(M, LX, ML, 2 * n, 13, 1);
	for (int i = 0; i < 2 * n; i++)
	{
		vv[i] = ML[i] - W[i];
	}
}

//求待定点的坐标仪坐标（改正系统误差）
void ImagePointCorrection(double &x_, double &y_, double x0, double y0, double *LX, PointData *Locate, int i)
{
	double r = 0.0;
	r = sqrt((Locate[i].x - x0)*(Locate[i].x - x0) + (Locate[i].y - y0)*(Locate[i].y - y0));//计算像点像径
	x_ = Locate[i].x + (Locate[i].x - x0)*r*r*LX[11] + (Locate[i].x - x0)*r*r*r*r*LX[12];
	y_ = Locate[i].y + (Locate[i].y - y0)*r*r*LX[11] + (Locate[i].y - y0)*r*r*r*r*LX[12];//求出经像点坐标的非线性改正后的像点坐标
}

//待定点物方空间坐标近似值的解算
void ObjectApproximation(double &X, double &Y, double &Z, double *LL, double *LR, double x_L, double y_L, double x_R, double y_R) {
	double *N = new double[4 * 3];
	double *S = new double[3 * 1];
	double *Q = new double[4 * 1];
	double *NT = new double[3 * 4];
	double *NTN = new double[3 * 3];
	double *NTQ = new double[3 * 1];
	//求N矩阵
	N[0] = LL[0] + x_L * LL[8];
	N[1] = LL[1] + x_L * LL[9];
	N[2] = LL[2] + x_L * LL[10];
	N[3] = LL[4] + x_L * LL[8];
	N[4] = LL[5] + x_L * LL[9];
	N[5] = LL[6] + x_L * LL[10];

	N[6] = LR[0] + x_R * LR[8];
	N[7] = LR[1] + x_R * LR[9];
	N[8] = LR[2] + x_R * LR[10];
	N[9] = LR[4] + x_R * LR[8];
	N[10] = LR[5] + x_R * LR[9];
	N[11] = LR[6] + x_R * LR[10];
	//求Q矩阵
	Q[0] = -(LL[3] + x_L);
	Q[1] = -(LL[7] + y_L);
	Q[2] = -(LR[3] + x_R);
	Q[3] = -(LR[7] + y_R);
	//法方程计算
	MatrixTranposition(N, NT, 4, 3);
	MatrixMultiplication(NT, N, NTN, 3, 4, 3);
	MatrixInversion(NTN, 3);
	MatrixMultiplication(NT, Q, NTQ, 3, 4, 1);
	MatrixMultiplication(NTN, NTQ, S, 3, 3, 1);

	X = S[0];
	Y = S[1];
	Z = S[2];//物方空间坐标近似值
}

//待定点物方空间坐标精确值的解算
void ObjectAccuration(double &VxL, double &VyL, double &VxR, double &VyR, double &X, double &Y, double &Z, double Xj, double Yj, double Zj, double *LL, double *LR, double xL, double yL, double xR, double yR, int i) {
	double AL = 0.0;
	double AR = 0.0;
	double *N = new double[4 * 3];
	double *S = new double[3 * 1];
	double *Q = new double[4 * 1];
	double *NT = new double[3 * 4];
	double *NTN = new double[3 * 3];
	double *NTQ = new double[3 * 1];
	double NS[4 * 1] = { 0 };
	double X_temp = 0, Y_temp = 0, Z_temp = 0;
	X = Xj; Y = Yj; Z = Zj;
	do
	{
		X_temp = X; Y_temp = Y; Z_temp = Z;

		//计算左右两片的A
		AL = LL[8] * X + LL[9] * Y + LL[10] * Z + 1;
		AR = LR[8] * X + LR[10] * Y + LR[10] * Z + 1;

		//计算N矩阵
		N[0] = -(LL[0] + xL * LL[8]) / AL;
		N[1] = -(LL[1] + xL * LL[9]) / AL;
		N[2] = -(LL[2] + xL * LL[10]) / AL;

		N[3] = -(LL[4] + yL * LL[8]) / AL;
		N[4] = -(LL[5] + yL * LL[9]) / AL;
		N[5] = -(LL[6] + yL * LL[10]) / AL;

		N[6] = -(LR[0] + xR * LR[8]) / AR;
		N[7] = -(LR[1] + xR * LR[9]) / AR;
		N[8] = -(LR[2] + xR * LR[10]) / AR;

		N[9] = -(LR[4] + yR * LR[8]) / AR;
		N[10] = -(LR[5] + yR * LR[9]) / AR;
		N[11] = -(LR[6] + yR * LR[10]) / AR;


		//计算Q矩阵
		Q[0] = (LL[3] + xL) / AL;
		Q[1] = (LL[7] + yL) / AL;

		Q[2] = (LR[3] + xR) / AR;
		Q[3] = (LR[7] + yR) / AR;

		//法方程的计算
		MatrixTranposition(N, NT, 4, 3);
		MatrixMultiplication(NT, N, NTN, 3, 4, 3);
		MatrixInversion(NTN, 3);
		MatrixMultiplication(NT, Q, NTQ, 3, 4, 1);
		MatrixMultiplication(NTN, NTQ, S, 3, 3, 1);

		X = S[0];
		Y = S[1];
		Z = S[2];
	} while (fabs(X - X_temp) > 0.0001 || fabs(Y - Y_temp > 0.0001) || fabs(Z - Z_temp) > 0.0001);//两次计算的物方空间坐标绝对值之差小于0.0001则停止迭代

	//计算误差方程式
	MatrixMultiplication(N, S, NS, 4, 3, 1);
	VxL = NS[0] - Q[0];
	VyL = NS[1] - Q[1];
	VxR = NS[2] - Q[2];
	VyR = NS[3] - Q[3];
}

//求外方位元素
void ExtraneousElement(double *&Angle, double *L, double x0, double y0, double f) {
	double t[9] = { 0 };
	double l[3] = { 0 };
	double X[3] = { 0 };
	double tT[3 * 3] = { 0 };
	double tTt[3 * 3] = { 0 };
	double tTl[3 * 1] = { 0 };
	double r3 = 0;
	double a3 = 0, b3 = 0, c3 = 0;
	double A = 0, B = 0, C = 0;
	double ds = 0, db = 0, b1 = 0, b2 = 0;

	l[0] = -L[3];
	l[1] = -L[7];
	l[2] = -1.0;

	t[0] = L[0];
	t[1] = L[1];
	t[2] = L[2];
	t[3] = L[4];
	t[4] = L[5];
	t[5] = L[6];
	t[6] = L[8];
	t[7] = L[9];
	t[8] = L[10];

	//求出像片的方向余弦
	r3 = sqrt(L[8] * L[8] + L[9] * L[9] + L[10] * L[10]);
	a3 = L[8] / r3;
	b3 = L[9] / r3;
	c3 = L[10] / r3;

	//求出比例尺不一误差ds和坐标轴不垂直性误差dβ
	A = (L[0] * L[0] + L[1] * L[1] + L[2] * L[2]) / (r3*r3) - x0 * x0;
	B = (L[4] * L[4] + L[5] * L[5] + L[6] * L[6]) / (r3*r3) - y0 * y0;
	C = (L[0] * L[4] + L[1] * L[5] + L[2] * L[6]) / (r3*r3) - x0 * y0;
	ds = sqrt(A / B) - 1;
	db = asin(sqrt((C*C) / (B*A)));

	//求出方向余弦b1,b2
	b1 = (L[1] + b3 * x0 + (y0*b3 + L[5])*(1 + ds)*sin(db)) / f;
	b2 = (y0*b3 + L[5])*(1 + ds)*cos(db) / f;

	//求出三个外方位角元素phi,omega,kappa
	Angle[3] = atan(-a3 / c3);
	Angle[4] = asin(-b3);
	Angle[5] = atan(b1 / b2);

	//计算法方程式
	MatrixTranposition(t, tT, 3, 3);
	MatrixMultiplication(tT, t, tTt, 3, 3, 3);
	MatrixInversion(tTt, 3);
	MatrixMultiplication(tT, l, tTl, 3, 3, 1);
	MatrixMultiplication(tTt, tTl, X, 3, 3, 1);

	//计算三个外方位直线元素Xs,Ys,Zs
	Angle[0] = X[0];
	Angle[1] = X[1];
	Angle[2] = X[2];
}

//精度评定
void AccuracyDescription(double &m0, double *&M, double *Vp, double *MTM, int n)
{
	double sum = 0;
	//计算单位全中误差
	for (int i = 0; i < 2 * n; i++)
	{
		sum += Vp[i] * Vp[i];
	}
	m0 = sqrt(sum / (2 * n - 13));

	//计算每个未知数的中误差
	for (int i = 0; i < 13; i++)
	{
		M[i] = sqrt(MTM[i + i * 13])*m0;
	}
}

int main()
{
	/*************************读取IMG_2319、IMG_2423、控制点坐标文件****************************/
	//读取IMG_2419文件
	ifstream ifsLeftImage("IMG_2419.txt");
	if (!ifsLeftImage) {
		printf("Open Left Image Data Error!");
		system("pause");
	}
	vector<double>LeftImageData;
	double l;
	while (ifsLeftImage >> l) {
		LeftImageData.push_back(l);
	}
	ifsLeftImage.close();

	//读取IMG_2423文件
	ifstream ifsRightImage("IMG_2423.txt");
	if (!ifsRightImage) {
		printf("Open Right Image Data Error!");
		system("pause");
	}
	vector<double>RightImageData;
	double r;
	while (ifsRightImage >> r) {
		RightImageData.push_back(r);
	}
	ifsRightImage.close();
	
	//读取控制点坐标文件
	ifstream ifsControlPoints("控制点坐标.txt");
	if (!ifsControlPoints) {
		printf("Open Control Points Data Error!");
		system("pause");
	}
	vector<double>ControlPointsData;
	double c;
	while (ifsControlPoints >> c) {
		ControlPointsData.push_back(c);
	}
	ifsControlPoints.close();
	int ObjectNumber = ControlPointsData[0];

	/***************************************坐标转换及对应********************************************/
	//左片像点和控制点物方坐标的对应
	int ergodicLeft = 0;
	PointData LeftPoints[54];
	for (int i = 1; i < ObjectNumber * 4 - 1; i += 4) {
		for (int j = 0; j < 78 * 3; j += 3) {
			if (ControlPointsData[i] == LeftImageData[j]) {
				LeftPoints[ergodicLeft].PointsNumber = ControlPointsData[i];
				LeftPoints[ergodicLeft].x = LeftImageData[j + 1] * pixsize;
				LeftPoints[ergodicLeft].y = LeftImageData[j + 2] * pixsize;
				LeftPoints[ergodicLeft].X = ControlPointsData[i + 2];
				LeftPoints[ergodicLeft].Y = ControlPointsData[i + 3];
				LeftPoints[ergodicLeft].Z = -ControlPointsData[i + 1];
				ergodicLeft++;
			}
		}
	}
	
	//右片像点和控制点物方坐标的对应
	int ergodicRight = 0;
	PointData RightPoints[59];
	for (int i = 1; i < ObjectNumber * 4 - 1; i += 4) {
		for (int j = 0; j < 83 * 3; j += 3) {
			if (ControlPointsData[i] == RightImageData[j]) {
				RightPoints[ergodicRight].PointsNumber = ControlPointsData[i];
				RightPoints[ergodicRight].x = RightImageData[j + 1] * pixsize;
				RightPoints[ergodicRight].y = RightImageData[j + 2] * pixsize;
				RightPoints[ergodicRight].X = ControlPointsData[i + 2];
				RightPoints[ergodicRight].Y = ControlPointsData[i + 3];
				RightPoints[ergodicRight].Z = -ControlPointsData[i + 1];
				ergodicRight++;
			}
		}
	}

	//左片像点和右片像点的对应
	int ergodicCommon = 0;
	PointData CommonPoints[29];
	for (int i = 0; i < ergodicLeft; i++) {
		for (int j = 0; j < ergodicRight; j++) {
			if (LeftPoints[i].PointsNumber == RightPoints[j].PointsNumber) {
				CommonPoints[ergodicCommon].PointsNumber = LeftPoints[i].PointsNumber;
				CommonPoints[ergodicCommon].X = LeftPoints[i].X;
				CommonPoints[ergodicCommon].Y = LeftPoints[i].Y;
				CommonPoints[ergodicCommon].Z = LeftPoints[i].Z;
				CommonPoints[ergodicCommon].xl = LeftPoints[i].x;
				CommonPoints[ergodicCommon].yl = LeftPoints[i].y;
				CommonPoints[ergodicCommon].xr = RightPoints[i].x;
				CommonPoints[ergodicCommon].yr = RightPoints[i].y;
				ergodicCommon++;
			}
		}
	}

	/*******************************存储左右相片同名点及对应物方控制点******************************************/
	PointData *NewLeftPoints = new PointData[ergodicCommon];
	PointData *NewRightPoints = new PointData[ergodicCommon];
	
	//将左片上的同名点像点坐标及对应物方控制点坐标存储
	for (int i = 0; i < ergodicCommon; i++) {
		NewLeftPoints[i].PointsNumber = CommonPoints[i].PointsNumber;
		NewLeftPoints[i].x = CommonPoints[i].xl;
		NewLeftPoints[i].y = CommonPoints[i].yl;
		NewLeftPoints[i].X = CommonPoints[i].X;
		NewLeftPoints[i].Y = CommonPoints[i].Y;
		NewLeftPoints[i].Z = CommonPoints[i].Z;
	}

	//将右片上的同名点像点坐标及对应物方控制点坐标存储
	for (int i = 0; i < ergodicCommon; i++) {
		NewRightPoints[i].PointsNumber = CommonPoints[i].PointsNumber;
		NewRightPoints[i].x = CommonPoints[i].xr;
		NewRightPoints[i].y = CommonPoints[i].yr;
		NewRightPoints[i].X = CommonPoints[i].X;
		NewRightPoints[i].Y = CommonPoints[i].Y;
		NewRightPoints[i].Z = CommonPoints[i].Z;
	}

	/******************************************设定各种初值********************************************/
	//左右像片L系数的近似值
	double *IL = new double[11 * 1];
	double *IR = new double[11 * 1];

	//左右像片内方位元素的近似值
	double x0L = 0.0;
	double y0L = 0.0;
	double x0R = 0.0;
	double y0R = 0.0;

	//L系数和k1、k2的精确值
	double *ILA = new double[13 * 1];
	double *IRA = new double[13 * 1];

	//经像点坐标的非线性改正后的像点坐标
	double dxL = 0.0;
	double dyL = 0.0;
	double dxR = 0.0;
	double dyR = 0.0;

	double *SA = new double[ergodicCommon];
	double *S = new double[ergodicCommon];
	double *V = new double[ergodicCommon];
	
	//左右像片的主距
	double fL = 0.0;
	double fR = 0.0;

	//内方位元素的精确值
	double x0LA = 0.0;
	double y0LA = 0.0;
	double fLA = 0.0;
	double x0RA = 0.0;
	double y0RA = 0.0;
	double fRA = 0.0;

	//外方位元素
	double *DL = new double[6 * 1];
	double *DR = new double[6 * 1];

	double *mpl = new double[13 * 1];
	double *mpr = new double[13 * 1];

	/******************************************DLT解算********************************************/
	//左右像片L系数的近似值
	LiApproximation(IL, LeftPoints, ergodicLeft);
	LiApproximation(IR, RightPoints, ergodicRight);

	double *VpL = new double[2 * ergodicLeft];
	double *VpR = new double[2 * ergodicRight];
	double *MTML = new double[13 * 13];
	double *MTMR = new double[13 * 13];

	//求出内方位元素的近似值
	InteriorElement(x0L, y0L, fL, IL);
	InteriorElement(x0R, y0R, fR, IR);

	//Li系数和内方位元素精确值求解
	LiAccuration(VpL, MTML, ILA, x0LA, y0LA, fLA, IL, LeftPoints, x0L, y0L, fL, ergodicLeft);
	LiAccuration(VpR, MTMR, IRA, x0RA, y0RA, fRA, IR, RightPoints, x0R, y0R, fR, ergodicRight);

	//待定点"坐标仪坐标"改正，待定点物方空间坐标近似值、精确值解算
	for (int i = 0; i < ergodicCommon; i++) {
		ImagePointCorrection(dxL, dyL, x0L, y0L, ILA, NewLeftPoints, i);
		ImagePointCorrection(dxR, dyR, x0R, y0R, IRA, NewRightPoints, i);
		ObjectApproximation(SA[3 * i + 0], SA[3 * i + 1], SA[3 * i + 2], ILA, IRA, dxL, dyL, dxR, dyR);
		ObjectAccuration(V[2 * i + 0], V[2 * i + 1], V[2 * i + 2], V[2 * i + 3], S[3 * i + 0], S[3 * i + 1], S[3 * i + 2], SA[3 * i + 0], SA[3 * i + 1], SA[3 * i + 2], ILA, IRA, dxL, dyL, dxR, dyR, i);
	}

	//求出外方位元素
	ExtraneousElement(DL, ILA, x0LA, y0LA, fLA);
	ExtraneousElement(DR, IRA, x0RA, y0RA, fRA);

	double m0L = 0.0;
	double m0R = 0.0;
	double *mL = new double[13];
	double *mR = new double[13];

	AccuracyDescription(m0L, mL, VpL, MTML, ergodicLeft);
	AccuracyDescription(m0R, mR, VpR, MTMR, ergodicRight);

	double *DX = new double[ergodicCommon];
	double *DY = new double[ergodicCommon];
	double *DZ = new double[ergodicCommon];
	double *DM = new double[ergodicCommon];

	for (int i = 0; i < ergodicCommon; i++) {
		DX[i] = CommonPoints[i].X - S[3 * i + 0];
		DY[i] = CommonPoints[i].Y - S[3 * i + 1];
		DZ[i] = CommonPoints[i].Z - S[3 * i + 2];
		DM[i] = ((V[2 * i] * V[2 * i] + V[2 * i + 1] * V[2 * i + 1] + V[2 * i + 2] * V[2 * i + 2] + V[2 * i + 3] * V[2 * i + 3]) / (4 - 3));
	}

	/******************************************输出结果********************************************/
	FILE* Result;
	Result = fopen("DLT解算结果.txt", "w");
	fprintf(Result, "*************************IMG_2419左影像片*************************\n\n");
	fprintf(Result, "外方位元素\n");
	fprintf(Result, "Xs=%.10fmm\tYs=%.10fmm\tZs=%.10f%mm\n%Phi=%.10frad\tOmega=%.10frad\tKappa=%.10frad\n\n", DL[0], DL[1], DL[2], DL[3], DL[4], DL[5]);
	fprintf(Result, "内方位元素\n");
	fprintf(Result, "x0=%.10fmm\ty0=%.10fmm\tf=%.10fmm\n\n", x0LA, y0LA, fLA);
	fprintf(Result, "像点坐标中误差\n");
	fprintf(Result, "m0l=%.10fmm\n\n\n", m0L);

	fprintf(Result, "L系数\n");
	for (int i = 0; i < 11; i++)
	{
		fprintf(Result, "L%d=%.10fmm\n", i + 1, ILA[i]);
	}
	fprintf(Result, "k1=%.10fmm\n", ILA[11]);
	fprintf(Result, "k2=%.10fmm\n", ILA[12]);

	fprintf(Result, "\n未知数中误差\n");
	for (int i = 0; i < 11; i++)
	{
		fprintf(Result, "L%d的中误差为\t%.10fmm\n", i + 1, mL[i]);
	}
	fprintf(Result, "k1的中误差为\t%.10fmm\n", mL[11]);
	fprintf(Result, "k2的中误差为\t%.10fmm\n\n", mL[12]);

	fprintf(Result, "像点观测值残差\n");
	fprintf(Result, "点号\t\tVx\t\tVy\n");
	for (int i = 0; i < ergodicLeft; i++)
	{
		fprintf(Result, "%d\t%.10fmm\t%.10fmm\n", LeftPoints[i].PointsNumber, VpL[i * 2], VpL[i * 2 + 1]);
	}

	fprintf(Result, "\n\n");
	fprintf(Result, "*************************IMG_2423右影像片*************************\n\n");
	fprintf(Result, "外方位元素\n");
	fprintf(Result, "Xs=%.10fmm\tYs=%.10fmm\tZs=%.10f%mm\n%Phi=%.10frad\tOmega=%.10frad\tKappa=%.10frad\n\n", DR[0], DR[1], DR[2], DR[3], DR[4], DR[5]);
	fprintf(Result, "内方位元素\n");
	fprintf(Result, "x0=%.10fmm\ty0=%.10fmm\tf=%.10fmm\n\n", x0RA, y0RA, fRA);
	fprintf(Result, "像点坐标中误差\n");
	fprintf(Result, "m0r=%.10fmm\n\n\n", m0R);

	fprintf(Result, "L系数\n");
	for (int i = 0; i < 11; i++)
	{
		fprintf(Result, "L%d=%.10fmm\n", i + 1, IRA[i]);
	}
	fprintf(Result, "k1=%.10fmm\n", IRA[11]);
	fprintf(Result, "k2=%.10fmm\n", IRA[12]);

	fprintf(Result, "\n未知数中误差\n");
	for (int i = 0; i < 11; i++)
	{
		fprintf(Result, "L%d的中误差为\t%.10fmm\n", i + 1, mR[i]);
	}
	fprintf(Result, "k1的中误差为\t%.10fmm\n", mR[11]);
	fprintf(Result, "k2的中误差为\t%.10fmm\n\n", mR[12]);

	fprintf(Result, "像点观测值残差\n");
	fprintf(Result, "点号\t\tVx\t\tVy\n");
	for (int i = 0; i < ergodicRight; i++)
	{
		fprintf(Result, "%d\t%.10fmm\t%.10fmm\n", RightPoints[i].PointsNumber, VpR[i * 2], VpR[i * 2 + 1]);
	}
	fprintf(Result, "\n\n");

	fprintf(Result, "将所有点作为检查点\n\n");
	fprintf(Result, "点号\t\tDx\t\t\tDy\t\t\tDz\t\t\tm0\n");
	for (int i = 0; i < ergodicCommon; i++)
	{
		fprintf(Result, "%d\t%.10fmm\t%.10fmm\t%.10fmm\t%.10fmm\n", CommonPoints[i].PointsNumber, DX[i], DY[i], DZ[i], DM[i]);
	}
	fclose(Result);
}

// 运行程序: Ctrl + F5 或调试 >“开始执行(不调试)”菜单
// 调试程序: F5 或调试 >“开始调试”菜单

// 入门提示: 
//   1. 使用解决方案资源管理器窗口添加/管理文件
//   2. 使用团队资源管理器窗口连接到源代码管理
//   3. 使用输出窗口查看生成输出和其他消息
//   4. 使用错误列表窗口查看错误
//   5. 转到“项目”>“添加新项”以创建新的代码文件，或转到“项目”>“添加现有项”以将现有代码文件添加到项目
//   6. 将来，若要再次打开此项目，请转到“文件”>“打开”>“项目”并选择 .sln 文件
