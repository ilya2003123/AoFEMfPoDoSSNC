#pragma once

#include "Problem.h"
#include <iosfwd>
#include <string>
#include <vector>

std::string formatDouble(double value);
std::string linearEquationString(double freeTerm, double xCoeff);

void writeVector(std::ostream& out, const std::string& name, const std::vector<double>& values);
void writeMatrix(std::ostream& out, const std::string& name, const Matrix& matrix);
void writeReport(std::ostream& out, const ProblemDefinition& problem, const SolveResult& result);
