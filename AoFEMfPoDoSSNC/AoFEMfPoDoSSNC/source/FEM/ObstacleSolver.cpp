#include "ObstacleSolver.h"

#include "../Parser/Parser.h"
#include "../coefficient/coefficient.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace
{
	double evaluateAt(const Expression& expression, double x)
	{
		return eval(expression, x);
	}
}

ObstacleSolveResult solveObstacleProblem(const ObstacleProblemDefinition& problem)
{
	if (problem.intervals < 2)
	{
		throw std::runtime_error("Obstacle problem needs at least 2 intervals");
	}
	if (problem.length <= 0.0)
	{
		throw std::runtime_error("String length must be positive");
	}
	if (problem.maxIterations <= 0)
	{
		throw std::runtime_error("Maximum iteration count must be positive");
	}
	if (problem.toleranceScale <= 0.0)
	{
		throw std::runtime_error("Tolerance scale must be positive");
	}

	Parser pParser(problem.p.c_str());
	Parser qParser(problem.q.c_str());
	Parser fParser(problem.f.c_str());
	Parser obstacleParser(problem.obstacle.c_str());
	Expression pExpr = pParser.parse();
	Expression qExpr = qParser.parse();
	Expression fExpr = fParser.parse();
	Expression obstacleExpr = obstacleParser.parse();

	const int intervals = problem.intervals;
	const int nodeCount = intervals + 1;
	const int unknownCount = intervals - 1;
	const double h = problem.length / intervals;

	ObstacleSolveResult result;
	result.matrix.assign(unknownCount, std::vector<double>(unknownCount, 0.0));
	result.rhs.assign(unknownCount, 0.0);
	result.obstacleValues.resize(nodeCount, 0.0);
	result.finalCoeff.resize(nodeCount, 0.0);
	result.obstacleAtLeft = evaluateAt(obstacleExpr, 0.0);
	result.obstacleAtRight = evaluateAt(obstacleExpr, problem.length);

	if (result.obstacleAtLeft > 1e-12 || result.obstacleAtRight > 1e-12)
	{
		throw std::runtime_error("Obstacle must not be above fixed boundary values u(0)=u(l)=0");
	}

	for (int node = 0; node < nodeCount; ++node)
	{
		const double x = node * h;
		result.obstacleValues[node] = evaluateAt(obstacleExpr, x);
	}

	for (int interval = 0; interval < intervals; ++interval)
	{
		const double left = interval * h;
		const double right = (interval + 1) * h;

		for (int localI = 0; localI < 2; ++localI)
		{
			const int globalI = interval + localI;
			if (globalI == 0 || globalI == intervals)
			{
				continue;
			}
			const int row = globalI - 1;

			result.rhs[row] += integrateGauss5(left, right, [&](double x)
			{
				const double shapeI = localI == 0 ? (right - x) / h : (x - left) / h;
				return (*fExpr.func)(x) * shapeI;
			});

			for (int localJ = 0; localJ < 2; ++localJ)
			{
				const int globalJ = interval + localJ;
				if (globalJ == 0 || globalJ == intervals)
				{
					continue;
				}
				const int column = globalJ - 1;
				const double dShapeI = localI == 0 ? -1.0 / h : 1.0 / h;
				const double dShapeJ = localJ == 0 ? -1.0 / h : 1.0 / h;

				result.matrix[row][column] += integrateGauss5(left, right, [&](double x)
				{
					const double shapeI = localI == 0 ? (right - x) / h : (x - left) / h;
					const double shapeJ = localJ == 0 ? (right - x) / h : (x - left) / h;
					return (*pExpr.func)(x) * dShapeI * dShapeJ + (*qExpr.func)(x) * shapeI * shapeJ;
				});
			}
		}
	}

	std::vector<double> unknowns(unknownCount, 0.0);
	for (int i = 0; i < unknownCount; ++i)
	{
		const int node = i + 1;
		unknowns[i] = std::max(0.0, result.obstacleValues[node]);
	}

	const double tolerance = problem.toleranceScale * h;
	for (int iteration = 1; iteration <= problem.maxIterations; ++iteration)
	{
		double maxDelta = 0.0;

		for (int i = 0; i < unknownCount; ++i)
		{
			if (std::abs(result.matrix[i][i]) < 1e-14)
			{
				throw std::runtime_error("Zero diagonal in projected Gauss-Seidel solver");
			}

			double sigma = 0.0;
			for (int j = 0; j < unknownCount; ++j)
			{
				if (j != i)
				{
					sigma += result.matrix[i][j] * unknowns[j];
				}
			}

			const double candidate = (result.rhs[i] - sigma) / result.matrix[i][i];
			const double projected = std::max(candidate, result.obstacleValues[i + 1]);
			maxDelta = std::max(maxDelta, std::abs(projected - unknowns[i]));
			unknowns[i] = projected;
		}

		result.iterations = iteration;
		result.maxDelta = maxDelta;
		if (maxDelta <= tolerance)
		{
			result.converged = true;
			break;
		}
	}

	for (int i = 0; i < unknownCount; ++i)
	{
		result.finalCoeff[i + 1] = unknowns[i];
	}

	result.solutionByInterval.resize(intervals, { 0.0, 0.0 });
	for (int interval = 0; interval < intervals; ++interval)
	{
		const double left = interval * h;
		const double leftValue = result.finalCoeff[interval];
		const double rightValue = result.finalCoeff[interval + 1];
		const double slope = (rightValue - leftValue) / h;
		result.solutionByInterval[interval][0] = leftValue - slope * left;
		result.solutionByInterval[interval][1] = slope;
	}

	return result;
}
