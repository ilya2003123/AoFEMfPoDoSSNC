#include "Equation.h"


Equation::Equation()
{
	//m_coeff = new double[m_index + 1];
}

Equation::Equation(double k, double b)
{
	//m_coeff = new double[m_index + 1];
	m_coeff[0] = k;
	m_coeff[1] = b;
	m_coeff[2] = 0;
}

Equation::Equation(int index)
{
	m_index = index;
	//m_coeff = new double[m_index + 1];
}

Equation::Equation(const Equation& other)
{
	m_index = other.m_index;
	//m_coeff = new double[other.m_index + 1];
	for (int i = 0; i < other.m_index + 1; i++)
		m_coeff[i] = other.m_coeff[i];
}

Equation::~Equation()
{
	//delete[] m_coeff;
}



Equation& Equation::operator=(const Equation& other)
{
	m_index = other.m_index;
	//m_coeff = new double[other.m_index + 1];
	for (int i = 0; i < other.m_index + 1; i++)
		m_coeff[i] = other.m_coeff[i];

	return *this;
}

std::string Equation::outputEquation() const
{
	std::string res;

	for (int i = 0; i < m_index + 1; i++)
	{
		if (m_coeff[i] != 0)
		{
			if (i >= 2)
			{
				res += std::to_string(static_cast<int>(round(m_coeff[i])));
				res += "x^" + std::to_string(i);
			}
			else if (i == 1)
			{
				if (m_coeff[i] > 0 && i != m_index)
				{
					res += "+";
				}
				res += std::to_string(static_cast<int>(round(m_coeff[i])));
				res += "x";
			}
			else if (i == 0)
			{
				if (m_coeff[i] > 0)
				{
					res += "+";
				}
				res += std::to_string(static_cast<int>(round(m_coeff[i])));
			}
		}
	}

	// Обработка случая, когда строка пуста
	if (res.empty())
	{
		res = "0";
	}

	return res;
}


Equation operator+(Equation& eq1, Equation& eq2)
{
	Equation eq3;
	for (int i = 0; i <= eq3.m_index; i++)
	{
		eq3.m_coeff[i] = 0;
	}

	for (int i = 0; i <= eq3.m_index; i++)
	{
		eq3.m_coeff[i] += eq1.m_coeff[i];
		eq3.m_coeff[i] += eq2.m_coeff[i];
	}

	return eq3;
}

Equation operator-(Equation& eq1, Equation& eq2)
{
	Equation eq3;
	for (int i = 0; i <= eq3.m_index; i++)
	{
		eq3.m_coeff[i] = 0;
	}

	for (int i = 0; i <= eq3.m_index; i++)
	{
		eq3.m_coeff[i] += eq1.m_coeff[i];
		eq3.m_coeff[i] -= eq2.m_coeff[i];
	}

	return eq3;
}

Equation operator*(Equation& eq1, Equation& eq2)
{
	const int p = eq1.m_index + eq2.m_index;
	Equation eq3(p);
	for (int i = 0; i <= p; i++)
	{
		eq3.m_coeff[i] = 0;
	}

	for (int i = 0; i <= eq1.m_index; i++)
	{
		for (int j = 0; j <= eq2.m_index; j++)
		{
			eq3.m_coeff[i + j] += eq1.m_coeff[i] * eq2.m_coeff[j];
		}
	}

	return eq3;
}
