#pragma once

#include "../FEM/Problem.h"
#include <string>

void printExpressionHelp();
ProblemDefinition problemFromConsole();
ProblemDefinition problemFromExactSolution();
int chooseRunMode();
