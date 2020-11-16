#include<iostream>
#include<cmath>
#include <Eigen/Eigen>
#include<Eigen/Core>
#include<Eigen/Dense>
using namespace std;
using namespace Eigen;
#define nR 3
#define nT 3
int main()
{
	int i, j, k;
	double temp;
	MatrixXcd H(nR, nT);
	typedef Matrix<complex<double>, nR, 1> VectornRcd;
	typedef Matrix<complex<double>, nT, 1> VectornTcd;
	VectornRcd x;
	VectornTcd c, y;
	HouseholderQR<MatrixXcd> qr;
	complex<double> t;
	MatrixXcd R, Q;
	cout << "Please input the real part of H:" << endl;
	for (i = 0; i < nR; i++)
	{
		for (j = 0; j < nT; j++)
		{
			cin >> temp;
			H.real()(i, j) = temp;
		}
	}
	cout << "Please input the imaginary part of H:" << endl;
	for (i = 0; i < nR; i++)
	{
		for (j = 0; j < nT; j++)
		{
			cin >> temp;
			H.imag()(i, j) = temp;
		}
	}
	cout << "H = " << endl << H << endl;
	cout << "Please input the real part of x:" << endl;
	for (i = 0; i < nR; i++)
	{
		cin >> temp;
		x.real()(i) = temp;
	}
	cout << "Please input the imaginary part of x:" << endl;
	for (i = 0; i < nR; i++)
	{
		cin >> temp;
		x.imag()(i) = temp;
	}
	qr.compute(H);
	R = qr.matrixQR().triangularView<Upper>();
	Q = qr.householderQ();//compute QR Decomposition
	y = Q.adjoint() * x;//y = Q^H*x
	c(nT - 1) = y(nT - 1) / R(nT - 1, nT - 1);
	if (abs(c.real()(nT - 1)) >= abs(c.imag()(nT - 1)))
	{
		c.real()(nT - 1) = 1;
		c.imag()(nT - 1) = 0;
	}
	else
	{
		c.real()(nT - 1) = 0;
		c.imag()(nT - 1) = 1;
	}//compute c_nT
	for (i = nT - 2; i >= 0; i--)
	{
		t = 0;
		for (j = i; j < nT; j++)
		{	
			t += R(i, j) * c(j);

		}
		c(i) = y(i) - t;
		c(i) /= R(i, i);
		if (abs(c.real()(i)) >= abs(c.imag()(i)))
		{
			c.real()(i) = 1;
			c.imag()(i) = 0;
		}
		else
		{
			c.real()(i) = 0;
			c.imag()(i) = 1;
		}
	}//solve linear equation
	cout << "c =" << endl << c << endl;
}