#pragma once

#include "Problem.h"
#include <array>
#include <string>

double evaluatePiece(const std::array<double, 2>& equation, double x);
void writePlotFile(const std::string& outputPath, const ProblemDefinition& problem, const SolveResult& result);
