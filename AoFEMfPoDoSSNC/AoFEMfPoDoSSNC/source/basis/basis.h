#pragma once

#include <iostream>
#include "Equation.h"

struct Point
{
	double x, y;
	friend std::ostream& operator<<(std::ostream& out, Point& point);
};



double linealCoeff(Point p1, Point p2);
double freeMember(Point p1, Point p2);

struct Cap
{
	Point m_leftBorder;
	Point m_rightBorder;
	Equation m_equation;

	Cap() = default;
	Cap(Point leftBorder, Point rightBorder)
		: m_leftBorder(leftBorder), m_rightBorder(rightBorder),
		m_equation(linealCoeff(leftBorder, rightBorder), freeMember(leftBorder, rightBorder))
	{ }

	friend std::ostream& operator<<(std::ostream& out, Cap& cap);
};

void line(Point p1, Point p2);
void equationStraightLine(Cap** omega, int m);
