#pragma once

#include "Problem.h"
#include <string>

SolveResult solveProblem(const ProblemDefinition& problem);

SolveResult calculate(const ProblemDefinition& problem, const std::string& outputPath = "output.txt");

void calculate(
	int intervals,
	double limiter,
	const std::string& p,
	const std::string& q,
	const std::string& f,
	double length = 1.0);

void calculate(int intervals, double limiter, double p, double q, double f);
