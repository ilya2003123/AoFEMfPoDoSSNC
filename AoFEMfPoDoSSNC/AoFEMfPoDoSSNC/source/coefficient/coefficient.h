#pragma once

#include "../basis/basis.h"
#include "../Functions/functions.h"
#include <cmath>
#include <stdexcept>
#include <vector>

using Matrix = std::vector<std::vector<double>>;

inline double capValue(const Cap& cap, double x)
{
	return cap.m_equation.m_coeff[1] * x + cap.m_equation.m_coeff[0];
}

inline double capDerivative(const Cap& cap)
{
	return cap.m_equation.m_coeff[1];
}

template <typename Integrand>
double integrateGauss5(double left, double right, Integrand integrand)
{
	static const double nodes[5] = {
		-0.9061798459386640,
		-0.5384693101056831,
		 0.0,
		 0.5384693101056831,
		 0.9061798459386640
	};
	static const double weights[5] = {
		0.2369268850561891,
		0.4786286704993665,
		0.5688888888888889,
		0.4786286704993665,
		0.2369268850561891
	};

	const double center = 0.5 * (left + right);
	const double halfLength = 0.5 * (right - left);
	double sum = 0.0;

	for (int i = 0; i < 5; ++i)
	{
		const double x = center + halfLength * nodes[i];
		sum += weights[i] * integrand(x);
	}

	return halfLength * sum;
}

inline double productPiecewiseDerivativeFunction(Cap* phi, Cap* ksi, int m)
{
	double sum = 0.0;
	for (int interval = 0; interval < m; ++interval)
	{
		const double left = phi[interval].m_leftBorder.x;
		const double right = phi[interval].m_rightBorder.x;
		sum += capDerivative(phi[interval]) * capDerivative(ksi[interval]) * (right - left);
	}
	return sum;
}

inline Matrix integrateProduct(Cap** phi, functions::Abstract& p, functions::Abstract& q, int m)
{
	Matrix integral(m, std::vector<double>(m, 0.0));

	for (int i = 0; i < m; ++i)
	{
		for (int j = 0; j < m; ++j)
		{
			for (int interval = 0; interval < m; ++interval)
			{
				const Cap& phiI = phi[i][interval];
				const Cap& phiJ = phi[j][interval];
				const double left = phiI.m_leftBorder.x;
				const double right = phiI.m_rightBorder.x;
				const double dPhiI = capDerivative(phiI);
				const double dPhiJ = capDerivative(phiJ);

				integral[i][j] += integrateGauss5(left, right, [&](double x)
				{
					return p(x) * dPhiI * dPhiJ + q(x) * capValue(phiI, x) * capValue(phiJ, x);
				});
			}
		}
	}

	return integral;
}

inline Matrix integrateProduct(Cap** phi, double p, double q, int m)
{
	functions::Const pConst(p);
	functions::Const qConst(q);
	return integrateProduct(phi, pConst, qConst, m);
}

inline std::vector<double> rightSystemCoefficientC(int m, functions::Abstract& f, Cap** phi)
{
	std::vector<double> rightVector(m, 0.0);

	for (int i = 0; i < m; ++i)
	{
		for (int interval = 0; interval < m; ++interval)
		{
			const Cap& phiI = phi[i][interval];
			const double left = phiI.m_leftBorder.x;
			const double right = phiI.m_rightBorder.x;

			rightVector[i] += integrateGauss5(left, right, [&](double x)
			{
				return f(x) * capValue(phiI, x);
			});
		}
	}

	return rightVector;
}

inline std::vector<double> rightSystemCoefficientC(int m, double f, Cap** phi)
{
	functions::Const fConst(f);
	return rightSystemCoefficientC(m, fConst, phi);
}

inline std::vector<double> rightSystemCoefficientD(int m)
{
	std::vector<double> rightVector(m, 0.0);
	if (m > 0)
	{
		rightVector[m - 1] = 1.0;
	}

	return rightVector;
}

inline std::vector<double> thomasAlgorithm(const Matrix& A, const std::vector<double>& d)
{
	const int n = static_cast<int>(d.size());
	if (n == 0)
	{
		return {};
	}
	if (n == 1)
	{
		if (std::abs(A[0][0]) < 1e-14)
		{
			throw std::runtime_error("Zero diagonal in 1x1 system");
		}
		return { d[0] / A[0][0] };
	}

	std::vector<double> cPrime(n, 0.0);
	std::vector<double> dPrime(n, 0.0);
	std::vector<double> x(n, 0.0);

	if (std::abs(A[0][0]) < 1e-14)
	{
		throw std::runtime_error("Zero first diagonal in Thomas algorithm");
	}

	cPrime[0] = A[0][1] / A[0][0];
	dPrime[0] = d[0] / A[0][0];

	for (int i = 1; i < n; ++i)
	{
		const double denominator = A[i][i] - A[i][i - 1] * cPrime[i - 1];
		if (std::abs(denominator) < 1e-14)
		{
			throw std::runtime_error("Zero denominator in Thomas algorithm");
		}
		if (i < n - 1)
		{
			cPrime[i] = A[i][i + 1] / denominator;
		}
		dPrime[i] = (d[i] - A[i][i - 1] * dPrime[i - 1]) / denominator;
	}

	x[n - 1] = dPrime[n - 1];
	for (int i = n - 2; i >= 0; --i)
	{
		x[i] = dPrime[i] - cPrime[i] * x[i + 1];
	}

	return x;
}
