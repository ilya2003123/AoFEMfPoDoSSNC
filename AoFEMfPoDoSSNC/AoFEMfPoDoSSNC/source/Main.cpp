#pragma warning(disable: 26451)

#include "calculate.h"
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

ProblemDefinition diplomaExample(int number)
{
	ProblemDefinition problem;
	problem.intervals = 5;

	switch (number)
	{
	case 1:
		problem.name = "Example 1";
		problem.length = 2.0;
		problem.limiter = 0.5;
		problem.p = "piecewise(1, x+1, x+2)";
		problem.q = "piecewise(1, 3*x^2, 2*x)";
		problem.f = "piecewise(1, -3*x^4+4*x+2, x^2-3*x-0.5)";
		problem.exactU = "piecewise(1, -x^2, 0.5*x-1.5)";
		problem.pointQ = { {1.0, 3.0} };
		problem.pointF = { {1.0, -8.5} };
		break;
	case 2:
		problem.name = "Example 2";
		problem.length = 2.0;
		problem.limiter = 0.5;
		problem.p = "piecewise(1, x+1, x+2)";
		problem.q = "piecewise(1, 3*x^2, 2*x)";
		problem.f = "piecewise(1, 3*x^7-25*x^4-20*x^3, 1.5, 5*x^2-3*x-2.5, -2*x^3+13*x+4)";
		problem.exactU = "piecewise(1, x^5, 1.5, 2.5*x-1.5, -x^2+4.5)";
		problem.pointQ = { {1.0, 3.0} };
		problem.pointF = { {1.0, 5.5}, {1.5, 19.25} };
		break;
	case 3:
		problem.name = "Example 3";
		problem.length = 1.0;
		problem.limiter = 0.05;
		problem.p = "1";
		problem.q = "1";
		problem.f = "piecewise(0.25, x, 0.5, x^2-1.8125, 0.75, 0.4375, -1.55*x+1.6)";
		problem.exactU = "piecewise(0.25, x, 0.5, x^2+3/16, 0.75, 7/16, -1.55*x+1.6)";
		problem.pointQ = { {0.25, 1.0}, {0.5, 2.0}, {0.75, 7.0} };
		problem.pointF = { {0.25, 0.75}, {0.5, 1.875}, {0.75, 4.6125} };
		break;
	case 4:
		problem.name = "Example 4";
		problem.length = 2.0;
		problem.limiter = 1.0;
		problem.p = "piecewise(1, x+1, x+2)";
		problem.q = "piecewise(1, 3*x^2, 2*x)";
		problem.f = "piecewise(1, 3*x^3-1, 3*x^3-12*x^2+5*x)";
		problem.exactU = "piecewise(1, x, 1.5*x^2-6*x+5.5)";
		problem.pointQ = { {1.0, 3.0} };
		problem.pointF = { {1.0, 14.0} };
		break;
	case 5:
		problem.name = "Example 5";
		problem.length = 2.0;
		problem.limiter = 2.0;
		problem.p = "1+x^2";
		problem.q = "piecewise(1, 1, 2*x)";
		problem.f = "piecewise(0.5, -2-5*x^2, 1, 0.25, 1.5, 4*x^2-7.5*x, 2*x^3-14*x^2+18*x-2)";
		problem.exactU = "piecewise(0.5, x^2, 1, 0.25, 1.5, 2*x-1.75, x^2-4*x+5)";
		problem.pointQ = { {1.0, 5.0} };
		problem.pointF = { {0.5, 1.25}, {1.0, -2.75}, {1.5, 9.75} };
		break;
	default:
		return diplomaExample(1);
	}

	return problem;
}

ProblemDefinition problemFromDiplomaExamples()
{
	std::cout << "Choose example for plotting:\n";
	std::cout << "  1 - example 1, l=2, m=0.5\n";
	std::cout << "  2 - example 2, l=2, m=0.5\n";
	std::cout << "  3 - example 3, l=1, m=0.05\n";
	std::cout << "  4 - example 4, l=2, m=1\n";
	std::cout << "  5 - example 5, l=2, m=2\n";
	std::cout << "Example: ";

	int exampleNumber = 1;
	std::cin >> exampleNumber;
	std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

	return diplomaExample(exampleNumber);
}

ProblemDefinition problemFromConsole()
{
	ProblemDefinition problem;

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

std::string buildGeneratedRightSideExpression(
	const std::string& p,
	const std::string& q,
	const std::string& exactU)
{
	return "-(dx((" + p + ")*dx(" + exactU + ")))+(" + q + ")*(" + exactU + ")";
}

ProblemDefinition problemFromExactSolution()
{
	ProblemDefinition problem;
	problem.generatedRightSide = true;
	problem.name = "Generated from exact solution";

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
	std::cout << "  1 - choose one of 5 diploma examples\n";
	std::cout << "  2 - enter p(x), q(x), f(x) manually\n";
	std::cout << "  3 - enter exact u(x), p(x), q(x); generate f(x)\n";
	std::cout << "Mode: ";

	int mode = 1;
	std::cin >> mode;
	std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

	return mode;
}

void runPlotScript()
{
	const int exitCode = std::system("python plot_solution.py");
	if (exitCode != 0)
	{
		std::cout << "Plot script failed. You can run: python plot_solution.py" << std::endl;
	}
}

int main()
{
	setlocale(LC_ALL, "rus");
	setlocale(LC_NUMERIC, "C");

	try
	{
		const int mode = chooseRunMode();
		ProblemDefinition problem;
		if (mode == 2)
		{
			problem = problemFromConsole();
		}
		else if (mode == 3)
		{
			problem = problemFromExactSolution();
		}
		else
		{
			problem = problemFromDiplomaExamples();
		}

		std::ofstream out("output.txt");
		if (!out)
		{
			throw std::runtime_error("Cannot open output.txt");
		}

		SolveResult result = solveProblem(problem, out);
		writePlotFile("output_for_plot.txt", problem, result);
		out.close();
		runPlotScript();

		std::cout << "Done. See output.txt, output_for_plot.txt and solution_plot.png" << std::endl;
		return 0;
	}
	catch (const std::exception& error)
	{
		std::cerr << "Error: " << error.what() << std::endl;
		return 1;
	}
}
