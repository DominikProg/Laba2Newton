#include <iostream>
#include <stdlib.h>
#include <iomanip>
#include <math.h>
using namespace std;
const int n = 2; // кол-во уравнений
double* Ffunc (double[]); //вектор невязки
double** Jfunc (double[]); //матрица Якоби
double* Gauss (double** , double []);
double* minusF (double[]);
double* iterX (double[], double[]);
double Delta1 (double F[]);
double Delta2(double newx[], double x[]);
int main()
{
	const int NIT = 100;
	double x[n] = { 1, -1 };
	double e1 = 1e-9;
	double e2 = 1e-9;
	int k = 0;
	cout << setw(5) << 'N' << setw(15) << "delta1" << setw(15) << "delta2" << endl;
	while (k < NIT)
	{
		double* F = new double[n];
		F = Ffunc(x);
		double** J = new double* [n];
		for (int i = 0; i < (n); i++)
			J[i] = new double(n);
		J = Jfunc(x);
		double* mF = new double[n];
			mF = minusF(F);
		double* deltax = new double[n];
		deltax=Gauss(J, mF);
		double* newx = new double[n];
		newx = iterX(x, deltax);
		double D1 = Delta1(F);
		double D2 = Delta2(newx, x);
		cout << setw(5) << k << setw(15) << D1 << setw(15) << D2 << endl;
		for (int i=0; i<n; i++)
		x[i] = newx[i];
		if (D1 <= e1 && D2 <= e2)
			break;
		k++;
		/*delete[] F;
		for (int i = 0; i < n; i++)
			delete[] J[i];
		delete[] mF;
		delete[] deltax;
		delete[] newx;*/
	}
	for (int i = 0; i < n; i++)
		cout << x[i] << ' ';
}
double* Ffunc(double x[])
{
	double* F = new double[n];
	F[0] = cos(0.4 * x[1] + x[0] * x[0]) + x[1] * x[1] + x[0] * x[0] - 1.6;
	F[1] = 1.5 * x[0] * x[0] - x[1] * x[1] / 0.36 - 1;
	return F;
}
double** Jfunc(double x[])
{
	double** J = new double* [n];
	for (int i = 0; i < n; i++)
		J[i] = new double[n];
	J[0][0] = -2 * sin(0.4 * x[1] + x[0] * x[0]) * x[0] + 2 * x[0];
	J[0][1] = -0.8 * sin(0.4 * x[1] + x[0] * x[0]) + 2 * x[1];
	J[1][0] = 3 * x[0];
	J[1][1] = -x[1] / 0.18;
	return J;
}
double* Gauss(double** matrix,  double svobod[])
{
	int masper[n];
	int k = 0;
	double* masx = new double[n];
	for (int i = 0; i < n; i++)
		masx[i] = 0;
	for (int i = 0; i < n; i++)
		masper[i] = i;
	for (int k = 0; k < n; k++)
	{
		double max = 0;
		int M = 0;
		int p = 0;
		for (int i = k; i < n; i++)
		{
			int l = masper[i];
			if (abs(matrix[l][k]) <= max) continue;
			else
			{
				M = l;
				p = i;
				max = abs(matrix[l][k]);
			}
		}
		masper[p] = masper[k];
		masper[k] = M;
		double Main = matrix[M][k];
		/*if (Main == 0)
		{
			cout << "error";
			exit;
		}*/
		for (int j = k; j < n; j++)
			matrix[M][j] /= Main;
		svobod[M] /= Main;
		for (int i = k + 1; i < n; i++)
		{
			int l = masper[i];
			svobod[l] -= matrix[l][k] * svobod[M];
			for (int j = k + 1; j < n; j++)
				matrix[l][j] -= matrix[l][k] * matrix[M][j];
			matrix[l][k] = 0;
		}
		if (k == n - 1)
		{
			int l = masper[k];
			svobod[l] /= matrix[l][k];
			masx[k] = svobod[l];   
		}
	}
	for (int k = n - 2; k >= 0; k--)
	{
		int l = masper[k];
		double sum = 0;
		for (int j = k + 1; j < n; j++)
			sum += matrix[l][j] * masx[j];
		masx[k] = (svobod[l] - sum) / matrix[l][k];
	}
	return masx;
}
double* minusF(double F[])
{
	double* mF = new double[n];
	for (int i = 0; i < n; i++)
		mF[i] = -F[i];
	return mF;
}
double* iterX(double x[], double deltax[])
{
	double* newx = new double[n];
	for (int i = 0; i < n; i++)
		newx[i] = x[i] + deltax[i];
	return newx;
}
double Delta1(double F[])
{
	double D1 = 0;
	for (int i = 0; i < n; i++)
		if (D1 < abs(F[i])) D1 = abs(F[i]);
	return D1;
}
double Delta2(double newx[], double x[])
{
	double D2 = 0;
	for (int i = 0; i < n; i++)
		if (abs(newx[i]) < 1)
		{
			if (D2 < abs(newx[i] - x[i])) D2 = abs(newx[i] - x[i]);
		}
		else
		{
			if (D2 < abs((newx[i] - x[i]) / newx[i])) D2 = abs((newx[i] - x[i]) / newx[i]);
		}
	return D2;
}
