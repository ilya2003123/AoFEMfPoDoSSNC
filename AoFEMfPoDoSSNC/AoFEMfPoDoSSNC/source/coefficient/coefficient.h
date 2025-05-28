#pragma once

#include "../basis/basis.h"
#include "cmath"

double productPiecewiseDerivativeFunction(Cap* phi, Cap* ksi, int m)
{
	double sum = 0.0;
	for (int i = 0; i < m; i++)
	{
		sum += phi[i].m_equation.m_coeff[1] * ksi[i].m_equation.m_coeff[1];
	}
	return sum;
}

std::vector<std::vector<double>> integrateProduct(Cap** phi, double p, double q, int m)
{
	// Интеграл от p*φ₁'*φ₂' + q*φ₁*φ₂ на [x1, x2]
	std::vector<std::vector<double>> integral;
	integral.resize(m, std::vector<double>(m, 0.0));


	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < m; j++)
		{
			double k1 = phi[i]->m_equation.m_coeff[1];
			double k2 = phi[j]->m_equation.m_coeff[1];

			double delta = (phi[i][j].m_rightBorder.x - phi[i][j].m_leftBorder.x);

			// Жёсткостная часть: ∫p*φ₁'*φ₂'dx = p*k1*k2*(x2 - x1)
			integral[i][j] += p * productPiecewiseDerivativeFunction(phi[i], phi[j], m) * delta;
		}
	}

	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < m; j++)
		{
			// Массовая часть: ∫q*φ₁*φ₂dx, φ₁(x) = a1*x + b1, φ₂(x) = a2*x + b2
			for (int k = 0; k < m + 1; k++)
			{
				double a1 = phi[i][k].m_equation.m_coeff[1];
				double b1 = phi[i][k].m_equation.m_coeff[0];
				double a2 = phi[j][k].m_equation.m_coeff[1];
				double b2 = phi[j][k].m_equation.m_coeff[0];

				double a11 = a1 * a2 / ((m + 1) * (m + 1) * (m + 1));
				double a12 = (a1 * b2 + a2 * b1) / ((m + 1) * (m + 1));
				double a13 = b1 * b2 * (phi[i][j].m_rightBorder.x - phi[i][j].m_leftBorder.x);

				integral[i][j] += q * (a11 / 3 + a12 / 2 + a13);
			}
		}
	}
	return integral;
}

std::vector<double> rightSystemCoefficientC(int m, double f, Cap** phi)
{
	std::vector<double> sum;
	sum.resize(m, 0.0);

	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < m + 1; j++)
		{
			double a = phi[i][j].m_equation.m_coeff[1];
			double b = phi[i][j].m_equation.m_coeff[0];

			double delta = phi[i][j].m_rightBorder.x - phi[i][j].m_leftBorder.x;

			sum[i] += (a * delta * delta) / 2 + b * delta;
		}
	}

	return sum;
}

std::vector<double> thomasAlgorithm(const std::vector<std::vector<double>>& A,
	const std::vector<double>& d) {
	int n = d.size();
	std::vector<double> c_prime(n, 0.0);
	std::vector<double> d_prime(n, 0.0);
	std::vector<double> x(n, 0.0);

	// Прямой ход
	c_prime[0] = A[0][1] / A[0][0];
	d_prime[0] = d[0] / A[0][0];

	for (int i = 1; i < n; ++i) {
		double denominator = A[i][i] - A[i][i - 1] * c_prime[i - 1];
		if (i < n - 1) {
			c_prime[i] = A[i][i + 1] / denominator;
		}
		d_prime[i] = (d[i] - A[i][i - 1] * d_prime[i - 1]) / denominator;
	}

	// Обратный ход
	x[n - 1] = d_prime[n - 1];
	for (int i = n - 2; i >= 0; --i) {
		x[i] = d_prime[i] - c_prime[i] * x[i + 1];
	}

	return x;
}

std::vector<double> rightSystemCoefficientD(int m)
{
	std::vector<double> rightVector;
	rightVector.resize(m - 1, 0.0);
	rightVector.push_back(1.0);

	return rightVector;
}