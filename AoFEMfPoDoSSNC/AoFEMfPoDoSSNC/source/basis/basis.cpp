#include "basis.h"

double linearCoefficient(Point p1, Point p2)
{
	double m = static_cast<double>((p2.y - p1.y) / (p2.x - p1.x));
	return m;
}

double constantTerm(Point p1, Point p2)
{
	double k = static_cast<double>((p2.y - p1.y) / (p2.x - p1.x));
	double b = static_cast<double>(p1.y - k * p1.x);
	return b;
}

void buildLinearBasis(Basis& phi, int m, double length)
{
	const double h = length / m;

	for (int i = 0; i < m; i++)
	{
		const int peakNode = i + 1;
		for (int j = 0; j < m; j++)
		{
			const double xLeft = j * h;
			const double xRight = (j + 1) * h;
			const double yLeft = (j == peakNode) ? 1.0 : 0.0;
			const double yRight = ((j + 1) == peakNode) ? 1.0 : 0.0;

			phi[i][j].m_leftBorder = { xLeft, yLeft };
			phi[i][j].m_rightBorder = { xRight, yRight };
		}
	}

	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < m; j++)
		{
			phi[i][j].m_equation =
				Equation(constantTerm(phi[i][j].m_leftBorder, phi[i][j].m_rightBorder),
					linearCoefficient(phi[i][j].m_leftBorder, phi[i][j].m_rightBorder));
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
