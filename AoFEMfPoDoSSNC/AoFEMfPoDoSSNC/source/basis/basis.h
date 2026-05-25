#pragma once

#include <iostream>
#include <vector>
#include "Equation.h"

struct Point
{
	double x, y;
	friend std::ostream& operator<<(std::ostream& out, Point& point);
};

double linearCoefficient(Point p1, Point p2);
double constantTerm(Point p1, Point p2);

struct Cap
{
	Point m_leftBorder;
	Point m_rightBorder;
	Equation m_equation;

	Cap() = default;
	Cap(Point leftBorder, Point rightBorder)
		: m_leftBorder(leftBorder), m_rightBorder(rightBorder),
		m_equation(constantTerm(leftBorder, rightBorder), linearCoefficient(leftBorder, rightBorder))
	{ }

	friend std::ostream& operator<<(std::ostream& out, Cap& cap);
};

using Basis = std::vector<std::vector<Cap>>;

void buildLinearBasis(Basis& phi, int m, double length = 1.0);
