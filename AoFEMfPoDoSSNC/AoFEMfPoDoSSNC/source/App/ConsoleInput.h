#pragma once

#include "../FEM/Problem.h"
#include <string>

void printExpressionHelp();
ProblemDefinition problemFromConsole();
ProblemDefinition problemFromExactSolution();
int chooseRunMode();
std::string buildGeneratedRightSideExpression(
	const std::string& p,
	const std::string& q,
	const std::string& exactU);

std::string buildGeneratedRightSideDisplayExpression(
	const std::string& p,
	const std::string& q,
	const std::string& exactU);
