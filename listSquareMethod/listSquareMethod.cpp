// listSquareMethod.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"

#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

//最小二乘拟合相关函数定义
double sum(vector<double> Vnum, int n);
double MutilSum(vector<double> Vx, vector<double> Vy, int n);
double RelatePow(vector<double> Vx, int n, int ex);
double RelateMutiXY(vector<double> Vx, vector<double> Vy, int n, int ex);
void EMatrix(vector<double> Vx, vector<double> Vy, int n, int ex, double coefficient[]);
void CalEquation(int exp, double coefficient[]);
double F(double c[], int l, int m);
double Em[6][4];

#define N 1e-13

int fitting()
{
	//double x[] = { 0.00,0.056,0.112,0.168,0.224,0.280,0.336,0.392,0.448,0.504,0.560,0.616,0.672,0.728,0.784,0.84,0.896,0.952,0.1008,0.1064,1.12 };
	//double y[] = { 0.00,1.66,3.31,4.96,6.6,8.22,9.82,11.4,12.94,14.43,15.86,17.22,18.5,19.69,20.79,21.79,22.7,23.53,24.25,24.87,25.4 };
	double x[] = { -4, 1, 3, 5, 7 };
	double y[] = { 24, 14, 52, 114, 200};
	double a, b, c, m1, m2, m3, z1, z2, z3; a = b = c = 0;
	double sumx = 0, sumx2 = 0, sumx3 = 0, sumx4 = 0, sumy = 0, sumxy = 0, sumx2y = 0;
	int size = sizeof(x) / sizeof(double);
	for (int i = 0; i<size; i++)
	{
		sumx += x[i]; sumy += y[i];
		sumx2 += pow(x[i], 2); sumxy += x[i] * y[i];
		sumx3 += pow(x[i], 3); sumx2y += pow(x[i], 2)*y[i];
		sumx4 += pow(x[i], 4);
	}
	do {
		m1 = a; a = (sumx2y - sumx3*b - sumx2*c) / sumx4; z1 = (a - m1)*(a - m1);
		m2 = b; b = (sumxy - sumx*c - sumx3*a) / sumx2; z2 = (b - m2)*(b - m2);
		m3 = c; c = (sumy - sumx2*a - sumx*b) / 10; z3 = (c - m3)*(c - m3);
	} while ((z1>N) || (z2>N) || (z3>N));
	printf("a=%9.6f,\nb=%9.6f,\nc=%9.6f\n", a, b, c);
	printf("拟合方程为   y=%9.6fx*x+%9.6fx+%9.6f\n", a, b, c);
	return 0;
}


//主函数，这里将数据拟合成二次曲线
int main(int argc, char* argv[]) 
{

	//fitting();
	//double arry1[5] = { 0,0.25,0,5,0.75 };
	//double arry2[5] = { 1,1.283,1.649,2.212,2.178 };
	//double arry1[] = { -4, 1, 3, 5, 7 };
	//double arry2[] = { 24, 14, 52, 114, 200 };
	
	double arry1[] = { 13.795852, 3.721394, 3.647097, 3.572431, 3.497866, 3.428838, 3.354549, 3.280093, 3.205521, 3.131145, 3.06186, 2.987387, 2.912888, 2.838291, 2.763912, 2.702192, 2.627968, 2.553367, 2.478922, 2.404432, 2.335287, 2.260955, 2.186402, 2.112035, 2.037594, 1.968498, 1.894094, 1.819674, 1.745217, 1.676018 };
	double arry2[] = { 0.075581, 0.077766, 0.07991, 0.081958, 0.083886, 0.08561, 0.08726, 0.088734, 0.090136, 0.091222, 0.092055, 0.092583, 0.092873, 0.092899, 0.092616, 0.092192, 0.091447, 0.09047, 0.089243, 0.087747, 0.086081, 0.084165, 0.082127, 0.080019, 0.077874, 0.075489, 0.073247, 0.070856, 0.068929, 0.067039 };
	int arrySize = sizeof(arry1) / sizeof(double);
	double coefficient[5];
	memset(coefficient, 0, sizeof(double) * 5);
	vector<double> vx, vy;
	for (int i = 0; i < arrySize; i++)
	{
		vx.push_back(arry1[i]);
		vy.push_back(arry2[i]);
	}
	EMatrix(vx, vy, arrySize, 3, coefficient);
	printf("拟合方程为：y = %lf + %lfx + %lfx^2 \n", coefficient[1], coefficient[2], coefficient[3]);
	return 0;
}
//累加
double sum(vector<double> Vnum, int n)
{
	double dsum = 0;
	for (int i = 0; i < n; i++)
	{
		dsum += Vnum[i];
	}
	return dsum;
}
//乘积和
double MutilSum(vector<double> Vx, vector<double> Vy, int n)
{
	double dMultiSum = 0;
	for (int i = 0; i < n; i++)
	{
		dMultiSum += Vx[i] * Vy[i];
	}
	return dMultiSum;
}
//ex次方和
double RelatePow(vector<double> Vx, int n, int ex)
{
	double ReSum = 0;
	for (int i = 0; i < n; i++)
	{
		ReSum += pow(Vx[i], ex);
	}
	return ReSum;
}
//x的ex次方与y的乘积的累加
double RelateMutiXY(vector<double> Vx, vector<double> Vy, int n, int ex)
{
	double dReMultiSum = 0;
	for (int i = 0; i < n; i++)
	{
		dReMultiSum += pow(Vx[i], ex)*Vy[i];
	}
	return dReMultiSum;
}
//计算方程组的增广矩阵
void EMatrix(vector<double> Vx, vector<double> Vy, int n, int ex, double coefficient[])
{
	for (int i = 1; i <= ex; i++)
	{
		for (int j = 1; j <= ex; j++)
		{
			Em[i][j] = RelatePow(Vx, n, i + j - 2);
		}
		Em[i][ex + 1] = RelateMutiXY(Vx, Vy, n, i - 1);
	}
	Em[1][1] = n;
	CalEquation(ex, coefficient);
}
//求解方程
void CalEquation(int exp, double coefficient[])
{
	for (int k = 1; k < exp; k++) //消元过程
	{
		for (int i = k + 1; i < exp + 1; i++)
		{
			double p1 = 0;

			if (Em[k][k] != 0)
				p1 = Em[i][k] / Em[k][k];

			for (int j = k; j < exp + 2; j++)
				Em[i][j] = Em[i][j] - Em[k][j] * p1;
		}
	}
	coefficient[exp] = Em[exp][exp + 1] / Em[exp][exp];
	for (int l = exp - 1; l >= 1; l--)   //回代求解
		coefficient[l] = (Em[l][exp + 1] - F(coefficient, l + 1, exp)) / Em[l][l];
}
//供CalEquation函数调用
double F(double c[], int l, int m)
{
	double sum = 0;
	for (int i = l; i <= m; i++)
		sum += Em[l - 1][i] * c[i];
	return sum;
}
