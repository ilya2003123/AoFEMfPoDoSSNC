#include "ConsoleInput.h"

#include "../FEM/GeneratedRightSide.h"

#include <iostream>
#include <limits>
#include <stdexcept>

void printExpressionHelp()
{
	std::cout << "Expression syntax: use x, +, -, *, /, ^ and functions like sin(...), cos(...), exp(...), log(...).\n";
	std::cout << "Piecewise syntax: pw(border1, expr1, border2, expr2, lastExpr).\n";
	std::cout << "The longer name piecewise(...) is also supported.\n";
	std::cout << "It means expr1 for x <= border1, expr2 for x <= border2, otherwise lastExpr.\n\n";
}

ProblemDefinition problemFromConsole()
{
	ProblemDefinition problem;
	printExpressionHelp();

	std::cout << "N: ";
	std::cin >> problem.intervals;

	std::cout << "l: ";
	std::cin >> problem.length;

	std::cout << "m: ";
	std::cin >> problem.limiter;
	std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

	std::cout << "p(x): ";
	std::getline(std::cin, problem.p);

	std::cout << "q(x): ";
	std::getline(std::cin, problem.q);

	std::cout << "f(x): ";
	std::getline(std::cin, problem.f);

	std::cout << "Exact u(x), empty if unknown: ";
	std::getline(std::cin, problem.exactU);

	return problem;
}

ProblemDefinition problemFromExactSolution()
{
	ProblemDefinition problem;
	problem.generatedRightSide = true;
	problem.name = "Generated from exact solution";
	printExpressionHelp();

	std::cout << "N: ";
	std::cin >> problem.intervals;

	std::cout << "l: ";
	std::cin >> problem.length;

	std::cout << "m: ";
	std::cin >> problem.limiter;
	std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

	std::cout << "p(x): ";
	std::getline(std::cin, problem.p);

	std::cout << "q(x): ";
	std::getline(std::cin, problem.q);

	std::cout << "Exact u(x): ";
	std::getline(std::cin, problem.exactU);

	problem.f = buildGeneratedRightSideExpression(problem.p, problem.q, problem.exactU);
	return problem;
}

int chooseRunMode()
{
	std::cout << "Choose run mode:\n";
	std::cout << "  1 - enter p(x), q(x), f(x) manually\n";
	std::cout << "  2 - enter exact u(x), p(x), q(x); generate f(x)\n";
	std::cout << "Mode: ";

	int mode = 1;
	std::cin >> mode;
	std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

	if (mode != 1 && mode != 2)
	{
		throw std::runtime_error("Unknown run mode. Use 1 or 2.");
	}

	return mode;
}

