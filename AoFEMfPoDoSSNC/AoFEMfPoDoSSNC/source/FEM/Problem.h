#pragma once

#include "Types.h"
#include <array>
#include <string>
#include <vector>

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

struct ObstacleProblemDefinition
{
	std::string name;
	int intervals = 20;
	double length = 1.0;
	std::string p = "1";
	std::string q = "0";
	std::string f = "2";
	std::string obstacle = "2*x*(1-x)-0.1";
	int maxIterations = 10000;
	double toleranceScale = 1e-3;
};

struct ObstacleSolveResult
{
	Matrix matrix;
	std::vector<double> rhs;
	std::vector<double> obstacleValues;
	std::vector<double> finalCoeff;
	std::vector<std::array<double, 2>> solutionByInterval; // {free term, x coefficient}
	double obstacleAtLeft = 0.0;
	double obstacleAtRight = 0.0;
	int iterations = 0;
	bool converged = false;
	double maxDelta = 0.0;
};
