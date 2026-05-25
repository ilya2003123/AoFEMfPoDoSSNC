#pragma warning(disable: 26451)

#include "App/ConsoleInput.h"
#include "FEM/PlotDataWriter.h"
#include "FEM/ReportWriter.h"
#include "FEM/Solver.h"
#include <clocale>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <stdexcept>

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
		ProblemDefinition problem = (mode == 1)
			? problemFromConsole()
			: problemFromExactSolution();

		std::ofstream out("output.txt");
		if (!out)
		{
			throw std::runtime_error("Cannot open output.txt");
		}

		SolveResult result = solveProblem(problem);
		writeReport(out, problem, result);
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
