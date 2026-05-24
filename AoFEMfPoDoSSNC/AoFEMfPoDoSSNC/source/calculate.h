#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "Functions/functions.h"
#include "Operations/operations.h"
#include "Utils/Utils.h"
#include "Derivative/derivative.h"
#include "Parser/Parser.h"
#include "basis/basis.h"
#include "coefficient/coefficient.h"

using namespace utils;

struct ProblemDefinition
{
	std::string name;
	int intervals = 10;
	double length = 1.0;
	double limiter = 0.5;
	std::string p = "1";
	std::string q = "0";
	std::string f = "1";
	std::string exactU;
	bool generatedRightSide = false;
	std::vector<std::array<double, 2>> pointQ; // {x, jump of Q}
	std::vector<std::array<double, 2>> pointF; // {x, jump of F}
};

struct SolveResult
{
	Matrix matrix;
	std::vector<double> rhsC;
	std::vector<double> rhsD;
	std::vector<double> coeffC;
	std::vector<double> coeffD;
	std::vector<double> finalCoeff;
	std::vector<std::array<double, 2>> solutionByInterval; // {free term, x coefficient}
	double upAtEnd = 0.0;
	double wpAtEnd = 0.0;
	double correction = 0.0;
	std::string boundaryCase;
};

inline std::string formatDouble(double value)
{
	if (std::abs(value) < 1e-12)
	{
		value = 0.0;
	}

	std::ostringstream out;
	out << std::fixed << std::setprecision(10) << value;
	return out.str();
}

inline std::string linearEquationString(double freeTerm, double xCoeff)
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

inline void writeVector(std::ostream& out, const std::string& name, const std::vector<double>& values)
{
	out << name << ":\n";
	for (int i = 0; i < static_cast<int>(values.size()); ++i)
	{
		out << "  [" << i + 1 << "] = " << formatDouble(values[i]) << "\n";
	}
	out << "\n";
}

inline void writeMatrix(std::ostream& out, const std::string& name, const Matrix& matrix)
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

inline void writePointTerms(
	std::ostream& out,
	const std::string& name,
	const std::vector<std::array<double, 2>>& terms)
{
	if (terms.empty())
	{
		return;
	}

	out << name << ":\n";
	for (const auto& term : terms)
	{
		out << "  x = " << formatDouble(term[0])
			<< ", value = " << formatDouble(term[1]) << "\n";
	}
	out << "\n";
}

inline void writeBasis(std::ostream& out, Cap** phi, int intervals)
{
	out << "Basis functions phi_i on mesh intervals:\n";
	for (int i = 0; i < intervals; ++i)
	{
		out << "phi_" << i + 1 << ":\n";
		for (int j = 0; j < intervals; ++j)
		{
			const Cap& cap = phi[i][j];
			const double freeTerm = cap.m_equation.m_coeff[0];
			const double xCoeff = cap.m_equation.m_coeff[1];
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

inline double evaluatePiece(const std::array<double, 2>& equation, double x)
{
	return equation[1] * x + equation[0];
}

inline double evaluateBasis(Cap* phi, int intervals, double x)
{
	double value = 0.0;
	for (int interval = 0; interval < intervals; ++interval)
	{
		const Cap& cap = phi[interval];
		const bool inside = x >= cap.m_leftBorder.x - 1e-12 && x <= cap.m_rightBorder.x + 1e-12;
		if (!inside)
		{
			continue;
		}

		const double current = cap.m_equation.m_coeff[1] * x + cap.m_equation.m_coeff[0];
		if (std::abs(current) > std::abs(value))
		{
			value = current;
		}
	}

	if (std::abs(value) < 1e-12)
	{
		return 0.0;
	}
	return value;
}

inline void addPointMassTerms(
	Matrix& matrix,
	Cap** phi,
	int intervals,
	const std::vector<std::array<double, 2>>& terms)
{
	if (terms.empty())
	{
		return;
	}

	std::vector<double> values(intervals, 0.0);
	for (const auto& term : terms)
	{
		for (int i = 0; i < intervals; ++i)
		{
			values[i] = evaluateBasis(phi[i], intervals, term[0]);
		}

		for (int i = 0; i < intervals; ++i)
		{
			for (int j = 0; j < intervals; ++j)
			{
				matrix[i][j] += term[1] * values[i] * values[j];
			}
		}
	}
}

inline void addPointLoadTerms(
	std::vector<double>& rhs,
	Cap** phi,
	int intervals,
	const std::vector<std::array<double, 2>>& terms)
{
	if (terms.empty())
	{
		return;
	}

	for (const auto& term : terms)
	{
		for (int i = 0; i < intervals; ++i)
		{
			rhs[i] += term[1] * evaluateBasis(phi[i], intervals, term[0]);
		}
	}
}

inline void writePlotFile(
	const std::string& outputPath,
	const ProblemDefinition& problem,
	const SolveResult& result)
{
	std::ofstream out(outputPath);
	if (!out)
	{
		throw std::runtime_error("Cannot open plot output file: " + outputPath);
	}

	out << std::fixed << std::setprecision(10);
	out << "# AoFEMfPoDoSSNC plot data v1\n";
	if (!problem.name.empty())
	{
		out << "name " << problem.name << "\n";
	}
	out << "N " << problem.intervals << "\n";
	out << "l " << problem.length << "\n";
	out << "m " << problem.limiter << "\n";
	out << "p " << problem.p << "\n";
	out << "q " << problem.q << "\n";
	if (!problem.generatedRightSide)
	{
		out << "f " << problem.f << "\n";
	}
	if (!problem.exactU.empty())
	{
		out << "exact_u " << problem.exactU << "\n";
	}
	out << "case " << result.boundaryCase << "\n";
	out << "alpha " << result.correction << "\n";
	out << "up_l " << result.upAtEnd << "\n";
	out << "wp_l " << result.wpAtEnd << "\n\n";

	if (!problem.pointQ.empty())
	{
		out << "[point_q]\n";
		out << "# x jump_Q\n";
		for (const auto& term : problem.pointQ)
		{
			out << term[0] << " " << term[1] << "\n";
		}
		out << "\n";
	}

	if (!problem.pointF.empty())
	{
		out << "[point_f]\n";
		out << "# x jump_F\n";
		for (const auto& term : problem.pointF)
		{
			out << term[0] << " " << term[1] << "\n";
		}
		out << "\n";
	}

	out << "[coefficients]\n";
	out << "# i a_i\n";
	for (int i = 0; i < static_cast<int>(result.finalCoeff.size()); ++i)
	{
		out << i + 1 << " " << result.finalCoeff[i] << "\n";
	}
	out << "\n";

	out << "[pieces]\n";
	out << "# left right free_term x_coeff\n";
	for (int interval = 0; interval < problem.intervals; ++interval)
	{
		const double left = problem.length * interval / problem.intervals;
		const double right = problem.length * (interval + 1) / problem.intervals;
		out << left << " "
			<< right << " "
			<< result.solutionByInterval[interval][0] << " "
			<< result.solutionByInterval[interval][1] << "\n";
	}
	out << "\n";

	out << "[nodes]\n";
	out << "# x u_h(x)\n";
	for (int node = 0; node <= problem.intervals; ++node)
	{
		const double x = problem.length * node / problem.intervals;
		const int interval = std::min(node, problem.intervals - 1);
		const double value = evaluatePiece(result.solutionByInterval[interval], x);
		out << x << " " << value << "\n";
	}
}

inline SolveResult solveProblem(const ProblemDefinition& problem, std::ostream& out)
{
	if (problem.intervals <= 0)
	{
		throw std::runtime_error("Number of intervals must be positive");
	}
	if (problem.length <= 0.0)
	{
		throw std::runtime_error("String length must be positive");
	}
	if (problem.limiter <= 0.0)
	{
		throw std::runtime_error("Limiter m must be positive");
	}

	Parser pParser(problem.p.c_str());
	Parser qParser(problem.q.c_str());
	Parser fParser(problem.f.c_str());
	Expression pExpr = pParser.parse();
	Expression qExpr = qParser.parse();
	Expression fExpr = fParser.parse();

	const int n = problem.intervals;
	std::vector<std::vector<Cap>> phiStorage(n, std::vector<Cap>(n));
	std::vector<Cap*> phi(n, nullptr);
	for (int i = 0; i < n; ++i)
	{
		phi[i] = phiStorage[i].data();
	}

	equationStraightLine(phi.data(), n, problem.length);

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
	writePointTerms(out, "Point jumps of Q", problem.pointQ);
	writePointTerms(out, "Point jumps of F", problem.pointF);

	out << "Mesh nodes:\n";
	for (int node = 0; node <= n; ++node)
	{
		out << "  x_" << node << " = " << formatDouble(problem.length * node / n) << "\n";
	}
	out << "\n";

	writeBasis(out, phi.data(), n);

	SolveResult result;
	result.matrix = integrateProduct(phi.data(), *pExpr.func, *qExpr.func, n);
	result.rhsC = rightSystemCoefficientC(n, *fExpr.func, phi.data());
	addPointMassTerms(result.matrix, phi.data(), n, problem.pointQ);
	addPointLoadTerms(result.rhsC, phi.data(), n, problem.pointF);
	result.rhsD = rightSystemCoefficientD(n);

	writeMatrix(out, "Gram matrix A, A[i][k] = <phi_i, phi_k>", result.matrix);
	writeVector(out, "Right side b for u_p, b[k] = integral f * phi_k dx", result.rhsC);
	result.coeffC = thomasAlgorithm(result.matrix, result.rhsC);
	writeVector(out, "Coefficients c_i of u_p", result.coeffC);

	writeVector(out, "Right side e for w_p, e[k] = phi_k(l)", result.rhsD);
	result.coeffD = thomasAlgorithm(result.matrix, result.rhsD);
	writeVector(out, "Coefficients d_i of w_p", result.coeffD);

	result.upAtEnd = result.coeffC.back();
	result.wpAtEnd = result.coeffD.back();
	if (std::abs(result.wpAtEnd) < 1e-14)
	{
		throw std::runtime_error("w_p(l) is zero, boundary correction is impossible");
	}

	if (result.upAtEnd >= problem.limiter)
	{
		result.boundaryCase = "upper boundary";
		result.correction = (problem.limiter - result.upAtEnd) / result.wpAtEnd;
	}
	else if (result.upAtEnd <= -problem.limiter)
	{
		result.boundaryCase = "lower boundary";
		result.correction = (-problem.limiter - result.upAtEnd) / result.wpAtEnd;
	}
	else
	{
		result.boundaryCase = "free end";
		result.correction = 0.0;
	}

	out << "Boundary check:\n";
	out << "  u_p(l) = c_N = " << formatDouble(result.upAtEnd) << "\n";
	out << "  w_p(l) = d_N = " << formatDouble(result.wpAtEnd) << "\n";
	out << "  case   = " << result.boundaryCase << "\n";
	out << "  alpha  = " << formatDouble(result.correction) << "\n\n";

	result.finalCoeff.resize(n, 0.0);
	for (int i = 0; i < n; ++i)
	{
		result.finalCoeff[i] = result.coeffC[i] + result.correction * result.coeffD[i];
	}
	writeVector(out, "Final coefficients a_i = c_i + alpha * d_i", result.finalCoeff);

	result.solutionByInterval.resize(n, { 0.0, 0.0 });
	for (int interval = 0; interval < n; ++interval)
	{
		for (int i = 0; i < n; ++i)
		{
			result.solutionByInterval[interval][0] += result.finalCoeff[i] * phi[i][interval].m_equation.m_coeff[0];
			result.solutionByInterval[interval][1] += result.finalCoeff[i] * phi[i][interval].m_equation.m_coeff[1];
		}
	}

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
		const double value = evaluatePiece(result.solutionByInterval[interval], x);
		out << "  u_h(x_" << node << ") = u_h(" << formatDouble(x) << ") = " << formatDouble(value) << "\n";
	}
	out << "\n";

	return result;
}

inline SolveResult calculate(const ProblemDefinition& problem, const std::string& outputPath = "output.txt")
{
	std::ofstream outFile(outputPath);
	if (!outFile)
	{
		throw std::runtime_error("Cannot open output file: " + outputPath);
	}

	SolveResult result = solveProblem(problem, outFile);
	writePlotFile("output_for_plot.txt", problem, result);
	return result;
}

inline void calculate(int n, double m, const std::string& p, const std::string& q, const std::string& f, double length = 1.0)
{
	ProblemDefinition problem;
	problem.intervals = n;
	problem.length = length;
	problem.limiter = m;
	problem.p = p;
	problem.q = q;
	problem.f = f;
	calculate(problem);
}

inline void calculate(int n, double m, double p, double q, double f)
{
	calculate(n, m, std::to_string(p), std::to_string(q), std::to_string(f), 1.0);
}
