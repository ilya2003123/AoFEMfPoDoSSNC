#pragma once
#pragma warning(disable: 26451)

#include <string>
#include <iostream>

class Equation
{
public:
	Equation();
	Equation(double k, double b);
	Equation(int index);
	Equation(const Equation& other);
	~Equation();

public:
	friend Equation operator+(Equation& eq1, Equation& eq2);
	friend Equation operator-(Equation& eq1, Equation& eq2);
	friend Equation operator*(Equation& eq1, Equation& eq2);
	Equation& operator=(const Equation& eq);

	friend std::ostream& operator<<(std::ostream& out, Equation& eq)
	{
		out << eq.outputEquation();

		return out;
	}
public:
	std::string outputEquation() const;
	int m_index = 2;
	double m_coeff[3] = { 0.0 };
};


