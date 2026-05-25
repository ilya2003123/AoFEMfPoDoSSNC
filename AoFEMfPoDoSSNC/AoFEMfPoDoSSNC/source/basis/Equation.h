#pragma once
#pragma warning(disable: 26451)

#include <array>
#include <iostream>
#include <string>

class Equation
{
public:
	static constexpr int maxDegree = 2;
	using Coefficients = std::array<double, maxDegree + 1>;

	Equation();
	Equation(double constantTerm, double linearCoefficient);
	Equation(int degree);
	Equation(const Equation& other) = default;
	~Equation() = default;

public:
	friend Equation operator+(const Equation& eq1, const Equation& eq2);
	friend Equation operator-(const Equation& eq1, const Equation& eq2);
	friend Equation operator*(const Equation& eq1, const Equation& eq2);
	Equation& operator=(const Equation& eq) = default;

	friend std::ostream& operator<<(std::ostream& out, const Equation& eq)
	{
		out << eq.outputEquation();

		return out;
	}

public:
	std::string outputEquation() const;
	int degree() const;
	double coefficient(int degree) const;
	void setCoefficient(int degree, double value);
	double constant() const;
	double slope() const;

private:
	static void checkDegree(int degree);

	int m_degree = 0;
	Coefficients m_coeff{};
};


