#include "PlotDataWriter.h"

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <stdexcept>

double evaluatePiece(const std::array<double, 2>& equation, double x)
{
	return equation[1] * x + equation[0];
}

void writePlotFile(const std::string& outputPath, const ProblemDefinition& problem, const SolveResult& result)
{
	std::ofstream out(outputPath);
	if (!out)
	{
		throw std::runtime_error("Cannot open plot output file: " + outputPath);
	}

	out << std::fixed << std::setprecision(10);
	out << "# AoFEMfPoDoSSNC plot data v1\n";
	if (!problem.name.empty())
	{
		out << "name " << problem.name << "\n";
	}
	out << "N " << problem.intervals << "\n";
	out << "l " << problem.length << "\n";
	out << "m " << problem.limiter << "\n";
	out << "p " << problem.p << "\n";
	out << "q " << problem.q << "\n";
	if (!problem.f.empty())
	{
		out << "f " << problem.f << "\n";
	}
	if (!problem.exactU.empty())
	{
		out << "exact_u " << problem.exactU << "\n";
	}
	out << "case " << result.boundaryCase << "\n";
	out << "alpha " << result.correction << "\n";
	out << "up_l " << result.upAtEnd << "\n";
	out << "wp_l " << result.wpAtEnd << "\n\n";

	out << "[coefficients]\n";
	out << "# i a_i\n";
	for (int i = 0; i < static_cast<int>(result.finalCoeff.size()); ++i)
	{
		out << i + 1 << " " << result.finalCoeff[i] << "\n";
	}
	out << "\n";

	out << "[pieces]\n";
	out << "# left right free_term x_coeff\n";
	for (int interval = 0; interval < problem.intervals; ++interval)
	{
		const double left = problem.length * interval / problem.intervals;
		const double right = problem.length * (interval + 1) / problem.intervals;
		out << left << " "
			<< right << " "
			<< result.solutionByInterval[interval][0] << " "
			<< result.solutionByInterval[interval][1] << "\n";
	}
	out << "\n";

	out << "[nodes]\n";
	out << "# x u_h(x)\n";
	for (int node = 0; node <= problem.intervals; ++node)
	{
		const double x = problem.length * node / problem.intervals;
		const int interval = std::min(node, problem.intervals - 1);
		const double value = evaluatePiece(result.solutionByInterval[interval], x);
		out << x << " " << value << "\n";
	}
}
