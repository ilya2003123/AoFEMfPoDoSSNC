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
