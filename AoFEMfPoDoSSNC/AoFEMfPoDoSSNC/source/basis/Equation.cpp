#include "Equation.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <stdexcept>

namespace
{
	std::string formatCoefficient(double value)
	{
		std::ostringstream out;
		out << std::setprecision(10) << value;
		return out.str();
	}

	int actualDegree(const Equation::Coefficients& coefficients)
	{
		for (int degree = Equation::maxDegree; degree > 0; --degree)
		{
			if (std::abs(coefficients[degree]) > 1e-14)
			{
				return degree;
			}
		}
		return 0;
	}
}

Equation::Equation()
{
}

Equation::Equation(double constantTerm, double linearCoefficient)
{
	m_coeff[0] = constantTerm;
	m_coeff[1] = linearCoefficient;
	m_degree = actualDegree(m_coeff);
}

Equation::Equation(int degree)
	: m_degree(degree)
{
	checkDegree(degree);
}

std::string Equation::outputEquation() const
{
	std::string result;

	for (int degree = 0; degree <= m_degree; ++degree)
	{
		const double value = m_coeff[degree];
		if (std::abs(value) <= 1e-14)
		{
			continue;
		}

		if (!result.empty() && value > 0.0)
		{
			result += "+";
		}

		result += formatCoefficient(value);
		if (degree == 1)
		{
			result += "x";
		}
		else if (degree >= 2)
		{
			result += "x^" + std::to_string(degree);
		}
	}

	if (result.empty())
	{
		result = "0";
	}

	return result;
}

int Equation::degree() const
{
	return m_degree;
}

double Equation::coefficient(int degree) const
{
	checkDegree(degree);
	return m_coeff[degree];
}

void Equation::setCoefficient(int degree, double value)
{
	checkDegree(degree);
	m_coeff[degree] = value;
	m_degree = actualDegree(m_coeff);
}

double Equation::constant() const
{
	return coefficient(0);
}

double Equation::slope() const
{
	return coefficient(1);
}

void Equation::checkDegree(int degree)
{
	if (degree < 0 || degree > maxDegree)
	{
		throw std::out_of_range("Equation degree is outside the supported range");
	}
}

Equation operator+(const Equation& eq1, const Equation& eq2)
{
	Equation result(std::max(eq1.degree(), eq2.degree()));
	for (int degree = 0; degree <= Equation::maxDegree; ++degree)
	{
		result.setCoefficient(degree, eq1.coefficient(degree) + eq2.coefficient(degree));
	}
	return result;
}

Equation operator-(const Equation& eq1, const Equation& eq2)
{
	Equation result(std::max(eq1.degree(), eq2.degree()));
	for (int degree = 0; degree <= Equation::maxDegree; ++degree)
	{
		result.setCoefficient(degree, eq1.coefficient(degree) - eq2.coefficient(degree));
	}
	return result;
}

Equation operator*(const Equation& eq1, const Equation& eq2)
{
	const int resultDegree = eq1.degree() + eq2.degree();
	Equation::checkDegree(resultDegree);

	Equation result(resultDegree);
	for (int i = 0; i <= eq1.degree(); ++i)
	{
		for (int j = 0; j <= eq2.degree(); ++j)
		{
			const int degree = i + j;
			result.setCoefficient(
				degree,
				result.coefficient(degree) + eq1.coefficient(i) * eq2.coefficient(j));
		}
	}
	return result;
}
