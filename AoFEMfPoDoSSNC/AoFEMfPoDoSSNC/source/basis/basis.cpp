#include "basis.h"
#include <iostream>

void line(Point p1, Point p2)
{
	double m = (p2.y - p1.y) / (p2.x - p1.x);
	double b = p1.y - m * p1.x;
	std::cout << "Equation (" << p1.x << ", " << p1.y << ") and (" << p2.x << ", " << p2.y <<
		") is y =" << m << "x + " << b << std::endl;
}

double linealCoeff(Point p1, Point p2)
{
	double m = static_cast<double>((p2.y - p1.y) / (p2.x - p1.x));
	return m;
}

double freeMember(Point p1, Point p2)
{
	double k = static_cast<double>((p2.y - p1.y) / (p2.x - p1.x));
	double b = static_cast<double>(p1.y - k * p1.x);
	return b;
}

void equationStraightLine(Cap** phi, int m)
{
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < m + 1; j++)
		{
				if (i == j)
				{
					phi[i][j].m_leftBorder = { j * (1.0 / (m + 1)), 0 };
					phi[i][j].m_rightBorder = { (j + 1) * (1.0 / (m + 1)), 1 };
				}
				else if (i == j - 1)
				{
					phi[i][j].m_leftBorder = { j * (1.0 / (m + 1)), 1 };
					phi[i][j].m_rightBorder = { (j + 1) * (1.0 / (m + 1)), 0 };
				}
				else
				{
					phi[i][j].m_leftBorder = { j * (1.0 / (m + 1)), 0 };
					phi[i][j].m_rightBorder = { (j + 1) * (1.0 / (m + 1)), 0 };
				}
		}
	}

	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < m + 1; j++)
		{
			if (j != m)
			{
				phi[i][j].m_equation =
					Equation(freeMember(phi[i][j].m_leftBorder, phi[i][j].m_rightBorder),
						linealCoeff(phi[i][j].m_leftBorder, phi[i][j].m_rightBorder));
				/*phi[i][j].m_polynom =
					Polynomial(freeMember(phi[i][j].m_leftBorder, phi[i][j].m_rightBorder),
						linealCoeff(phi[i][j].m_leftBorder, phi[i][j].m_rightBorder));*/
			}
		}
	}
}

std::ostream& operator<<(std::ostream& out, Point& point)
{
	out << "[" << point.x << ", " << point.y << "]";
	return out;
}

std::ostream& operator<<(std::ostream& out, Cap& cap)
{
	out << "Equation:" << cap.m_equation << "\n";

	return out;
}