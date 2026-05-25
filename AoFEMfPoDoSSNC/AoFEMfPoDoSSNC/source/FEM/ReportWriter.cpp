#include "ReportWriter.h"

#include "../basis/basis.h"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <ostream>
#include <sstream>
#include <vector>

std::string formatDouble(double value)
{
	if (std::abs(value) < 1e-12)
	{
		value = 0.0;
	}

	std::ostringstream out;
	out << std::fixed << std::setprecision(10) << value;
	return out.str();
}

std::string linearEquationString(double freeTerm, double xCoeff)
{
	std::ostringstream out;
	out << formatDouble(xCoeff) << " * x";
	if (freeTerm >= 0.0)
	{
		out << " + " << formatDouble(freeTerm);
	}
	else
	{
		out << " - " << formatDouble(-freeTerm);
	}
	return out.str();
}

void writeVector(std::ostream& out, const std::string& name, const std::vector<double>& values)
{
	out << name << ":\n";
	for (int i = 0; i < static_cast<int>(values.size()); ++i)
	{
		out << "  [" << i + 1 << "] = " << formatDouble(values[i]) << "\n";
	}
	out << "\n";
}

void writeMatrix(std::ostream& out, const std::string& name, const Matrix& matrix)
{
	out << name << ":\n";
	for (const auto& row : matrix)
	{
		out << "  ";
		for (double value : row)
		{
			out << std::setw(16) << formatDouble(value);
		}
		out << "\n";
	}
	out << "\n";
}

namespace
{
	void writeBasis(std::ostream& out, const Basis& phi, int intervals)
	{
		out << "Basis functions phi_i on mesh intervals:\n";
		for (int i = 0; i < intervals; ++i)
		{
			out << "phi_" << i + 1 << ":\n";
			for (int j = 0; j < intervals; ++j)
			{
				const Cap& cap = phi[i][j];
				const double freeTerm = cap.m_equation.constant();
				const double xCoeff = cap.m_equation.slope();
				const bool nonZero = std::abs(freeTerm) > 1e-12 || std::abs(xCoeff) > 1e-12;

				if (nonZero)
				{
					out << "  [" << formatDouble(cap.m_leftBorder.x) << ", "
						<< formatDouble(cap.m_rightBorder.x) << "]: "
						<< linearEquationString(freeTerm, xCoeff) << "\n";
				}
			}
			out << "\n";
		}
	}
}

void writeReport(std::ostream& out, const ProblemDefinition& problem, const SolveResult& result)
{
	const int n = problem.intervals;
	Basis phi(n, std::vector<Cap>(n));
	buildLinearBasis(phi, n, problem.length);

	out << "============================================================\n";
	out << "Problem\n";
	if (!problem.name.empty())
	{
		out << "name        = " << problem.name << "\n";
	}
	out << "intervals N = " << problem.intervals << "\n";
	out << "length l    = " << formatDouble(problem.length) << "\n";
	out << "limiter m   = " << formatDouble(problem.limiter) << "\n";
	out << "p(x)        = " << problem.p << "\n";
	out << "q(x)        = " << problem.q << "\n";
	if (problem.generatedRightSide)
	{
		out << "right side  = generated internally from exact u, p and q\n";
	}
	else
	{
		out << "f(x)        = " << problem.f << "\n";
	}
	if (!problem.exactU.empty())
	{
		out << "exact u(x)  = " << problem.exactU << "\n";
	}
	out << "\n";

	out << "Mesh nodes:\n";
	for (int node = 0; node <= n; ++node)
	{
		out << "  x_" << node << " = " << formatDouble(problem.length * node / n) << "\n";
	}
	out << "\n";

	writeBasis(out, phi, n);

	writeMatrix(out, "Gram matrix A, A[i][k] = <phi_i, phi_k>", result.matrix);
	writeVector(out, "Right side b for u_p, b[k] = integral f * phi_k dx", result.rhsC);
	writeVector(out, "Coefficients c_i of u_p", result.coeffC);
	writeVector(out, "Right side e for w_p, e[k] = phi_k(l)", result.rhsD);
	writeVector(out, "Coefficients d_i of w_p", result.coeffD);

	out << "Boundary check:\n";
	out << "  u_p(l) = c_N = " << formatDouble(result.upAtEnd) << "\n";
	out << "  w_p(l) = d_N = " << formatDouble(result.wpAtEnd) << "\n";
	out << "  case   = " << result.boundaryCase << "\n";
	out << "  alpha  = " << formatDouble(result.correction) << "\n\n";

	writeVector(out, "Final coefficients a_i = c_i + alpha * d_i", result.finalCoeff);

	out << "Piecewise solution u_h(x):\n";
	for (int interval = 0; interval < n; ++interval)
	{
		const double left = phi[0][interval].m_leftBorder.x;
		const double right = phi[0][interval].m_rightBorder.x;
		out << "  [" << formatDouble(left) << ", " << formatDouble(right) << "]: "
			<< linearEquationString(result.solutionByInterval[interval][0], result.solutionByInterval[interval][1]) << "\n";
	}
	out << "\n";

	out << "Nodal values of u_h:\n";
	for (int node = 0; node <= n; ++node)
	{
		const double x = problem.length * node / n;
		const int interval = std::min(node, n - 1);
		const double value =
			result.solutionByInterval[interval][1] * x + result.solutionByInterval[interval][0];
		out << "  u_h(x_" << node << ") = u_h(" << formatDouble(x) << ") = " << formatDouble(value) << "\n";
	}
	out << "\n";
}
