#pragma once

#include "../FEM/Problem.h"

#include <string>
#include <vector>

struct PlotSeries
{
	std::vector<double> femX;
	std::vector<double> femY;
	std::vector<double> nodeX;
	std::vector<double> nodeY;
	std::vector<double> exactX;
	std::vector<double> exactY;
	std::vector<double> obstacleX;
	std::vector<double> obstacleY;
	bool hasExact = false;
	bool hasObstacle = false;
	bool showLimiter = true;
};

struct SolutionPiece
{
	double left = 0.0;
	double right = 0.0;
	double freeTerm = 0.0;
	double xCoeff = 0.0;
};

struct LoadedPlot
{
	bool loaded = false;
	bool resetPlotView = true;
	std::string path;
	std::string boundaryCase;
	ProblemDefinition problem;
	PlotSeries plot;
};

struct GuiState
{
	int problemType = 1;
	int mode = 1;
	int intervals = 20;
	double length = 1.0;
	double limiter = 1.0;
	std::string p = "1";
	std::string q = "0";
	std::string f = "pw(0.5, -2, 2)";
	std::string exactU = "pw(0.5, x^2, -x^2+2*x-0.5)";
	bool updatePythonPlot = true;
	bool hasResult = false;
	bool resetPlotView = true;
	std::string status = "Ready";
	std::string generatedF;
	ProblemDefinition problem;
	SolveResult result;
	ObstacleProblemDefinition obstacleProblem;
	ObstacleSolveResult obstacleResult;
	std::string obstacleF = "10";
	std::string obstaclePsi = "2*x*(1-x)-0.1";
	int obstacleMaxIterations = 10000;
	double obstacleToleranceScale = 1e-3;
	PlotSeries plot;
	LoadedPlot comparisonPlots[2];
	float comparisonSplit = 0.5f;
	std::string comparisonStatus = "Load up to two saved plots.";
};

void initializeGuiState(GuiState& state);
void drawFemGui(GuiState& state);
