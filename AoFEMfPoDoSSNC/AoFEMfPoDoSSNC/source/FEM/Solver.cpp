#include "Solver.h"

#include "PlotDataWriter.h"
#include "ReportWriter.h"
#include "../Parser/Parser.h"
#include "../basis/basis.h"
#include "../coefficient/coefficient.h"
#include <cmath>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

SolveResult solveProblem(const ProblemDefinition& problem)
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
	Basis phi(n, std::vector<Cap>(n));
	buildLinearBasis(phi, n, problem.length);

	SolveResult result;
	result.matrix = assembleSystemMatrix(phi, *pExpr.func, *qExpr.func, n);
	result.rhsC = assembleLoadVector(n, *fExpr.func, phi);
	result.rhsD = assembleBoundaryVector(n);

	result.coeffC = solveTridiagonalSystem(result.matrix, result.rhsC);
	result.coeffD = solveTridiagonalSystem(result.matrix, result.rhsD);

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

	result.finalCoeff.resize(n, 0.0);
	for (int i = 0; i < n; ++i)
	{
		result.finalCoeff[i] = result.coeffC[i] + result.correction * result.coeffD[i];
	}

	result.solutionByInterval.resize(n, { 0.0, 0.0 });
	for (int interval = 0; interval < n; ++interval)
	{
		for (int i = 0; i < n; ++i)
		{
			result.solutionByInterval[interval][0] += result.finalCoeff[i] * phi[i][interval].m_equation.constant();
			result.solutionByInterval[interval][1] += result.finalCoeff[i] * phi[i][interval].m_equation.slope();
		}
	}

	return result;
}

SolveResult calculate(const ProblemDefinition& problem, const std::string& outputPath)
{
	std::ofstream outFile(outputPath);
	if (!outFile)
	{
		throw std::runtime_error("Cannot open output file: " + outputPath);
	}

	SolveResult result = solveProblem(problem);
	writeReport(outFile, problem, result);
	writePlotFile("output_for_plot.txt", problem, result);
	return result;
}

void calculate(
	int intervals,
	double limiter,
	const std::string& p,
	const std::string& q,
	const std::string& f,
	double length)
{
	ProblemDefinition problem;
	problem.intervals = intervals;
	problem.length = length;
	problem.limiter = limiter;
	problem.p = p;
	problem.q = q;
	problem.f = f;
	calculate(problem);
}

void calculate(int intervals, double limiter, double p, double q, double f)
{
	calculate(intervals, limiter, std::to_string(p), std::to_string(q), std::to_string(f), 1.0);
}
