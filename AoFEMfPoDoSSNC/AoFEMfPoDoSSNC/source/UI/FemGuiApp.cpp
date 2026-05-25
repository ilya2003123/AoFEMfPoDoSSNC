#define NOMINMAX

#include "FemGuiApp.h"

#include "../App/ConsoleInput.h"
#include "../FEM/PlotDataWriter.h"
#include "../FEM/ReportWriter.h"
#include "../FEM/Solver.h"
#include "../Parser/Parser.h"
#include "imgui.h"
#include "imgui_impl_dx11.h"
#include "imgui_impl_win32.h"
#include "imgui_stdlib.h"
#include "implot.h"

#include <d3d11.h>
#include <algorithm>
#include <commdlg.h>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <windows.h>

extern IMGUI_IMPL_API LRESULT ImGui_ImplWin32_WndProcHandler(HWND hwnd, UINT msg, WPARAM wParam, LPARAM lParam);

namespace
{
	ID3D11Device* gDevice = nullptr;
	ID3D11DeviceContext* gDeviceContext = nullptr;
	IDXGISwapChain* gSwapChain = nullptr;
	ID3D11RenderTargetView* gMainRenderTargetView = nullptr;

	struct PlotSeries
	{
		std::vector<double> femX;
		std::vector<double> femY;
		std::vector<double> nodeX;
		std::vector<double> nodeY;
		std::vector<double> exactX;
		std::vector<double> exactY;
		bool hasExact = false;
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
		PlotSeries plot;
		LoadedPlot comparisonPlots[2];
		float comparisonSplit = 0.5f;
		std::string comparisonStatus = "Load up to two saved plots.";
	};

	void createRenderTarget()
	{
		ID3D11Texture2D* backBuffer = nullptr;
		gSwapChain->GetBuffer(0, IID_PPV_ARGS(&backBuffer));
		gDevice->CreateRenderTargetView(backBuffer, nullptr, &gMainRenderTargetView);
		backBuffer->Release();
	}

	void cleanupRenderTarget()
	{
		if (gMainRenderTargetView)
		{
			gMainRenderTargetView->Release();
			gMainRenderTargetView = nullptr;
		}
	}

	bool createDeviceD3D(HWND hwnd)
	{
		DXGI_SWAP_CHAIN_DESC desc = {};
		desc.BufferCount = 2;
		desc.BufferDesc.Width = 0;
		desc.BufferDesc.Height = 0;
		desc.BufferDesc.Format = DXGI_FORMAT_R8G8B8A8_UNORM;
		desc.BufferDesc.RefreshRate.Numerator = 60;
		desc.BufferDesc.RefreshRate.Denominator = 1;
		desc.Flags = DXGI_SWAP_CHAIN_FLAG_ALLOW_MODE_SWITCH;
		desc.BufferUsage = DXGI_USAGE_RENDER_TARGET_OUTPUT;
		desc.OutputWindow = hwnd;
		desc.SampleDesc.Count = 1;
		desc.SampleDesc.Quality = 0;
		desc.Windowed = TRUE;
		desc.SwapEffect = DXGI_SWAP_EFFECT_DISCARD;

		const D3D_FEATURE_LEVEL featureLevelArray[2] = {
			D3D_FEATURE_LEVEL_11_0,
			D3D_FEATURE_LEVEL_10_0,
		};
		D3D_FEATURE_LEVEL featureLevel = {};

		const HRESULT result = D3D11CreateDeviceAndSwapChain(
			nullptr,
			D3D_DRIVER_TYPE_HARDWARE,
			nullptr,
			0,
			featureLevelArray,
			2,
			D3D11_SDK_VERSION,
			&desc,
			&gSwapChain,
			&gDevice,
			&featureLevel,
			&gDeviceContext);

		if (result != S_OK)
		{
			return false;
		}

		createRenderTarget();
		return true;
	}

	void cleanupDeviceD3D()
	{
		cleanupRenderTarget();
		if (gSwapChain)
		{
			gSwapChain->Release();
			gSwapChain = nullptr;
		}
		if (gDeviceContext)
		{
			gDeviceContext->Release();
			gDeviceContext = nullptr;
		}
		if (gDevice)
		{
			gDevice->Release();
			gDevice = nullptr;
		}
	}

	std::string trim(const std::string& value)
	{
		const std::size_t first = value.find_first_not_of(" \t\r\n");
		if (first == std::string::npos)
		{
			return "";
		}
		const std::size_t last = value.find_last_not_of(" \t\r\n");
		return value.substr(first, last - first + 1);
	}

	std::string wideToUtf8(const wchar_t* value)
	{
		const int needed = WideCharToMultiByte(CP_UTF8, 0, value, -1, nullptr, 0, nullptr, nullptr);
		if (needed <= 1)
		{
			return "";
		}

		std::vector<char> buffer(static_cast<std::size_t>(needed), '\0');
		WideCharToMultiByte(CP_UTF8, 0, value, -1, buffer.data(), needed, nullptr, nullptr);
		return std::string(buffer.data());
	}

	std::string showOpenPlotDialog()
	{
		wchar_t fileName[MAX_PATH] = L"";
		OPENFILENAMEW dialog = {};
		dialog.lStructSize = sizeof(dialog);
		dialog.lpstrFile = fileName;
		dialog.nMaxFile = MAX_PATH;
		dialog.lpstrFilter = L"FEM plot data (*.txt)\0*.txt\0All files (*.*)\0*.*\0";
		dialog.nFilterIndex = 1;
		dialog.Flags = OFN_PATHMUSTEXIST | OFN_FILEMUSTEXIST | OFN_NOCHANGEDIR;

		if (!GetOpenFileNameW(&dialog))
		{
			return "";
		}
		return wideToUtf8(fileName);
	}

	std::string showSavePlotDialog()
	{
		wchar_t fileName[MAX_PATH] = L"fem_plot.txt";
		OPENFILENAMEW dialog = {};
		dialog.lStructSize = sizeof(dialog);
		dialog.lpstrFile = fileName;
		dialog.nMaxFile = MAX_PATH;
		dialog.lpstrFilter = L"FEM plot data (*.txt)\0*.txt\0All files (*.*)\0*.*\0";
		dialog.lpstrDefExt = L"txt";
		dialog.nFilterIndex = 1;
		dialog.Flags = OFN_PATHMUSTEXIST | OFN_OVERWRITEPROMPT | OFN_NOCHANGEDIR;

		if (!GetSaveFileNameW(&dialog))
		{
			return "";
		}
		return wideToUtf8(fileName);
	}

	void appendExactSamples(PlotSeries& plot, const std::string& exactU, double length)
	{
		const int exactSamples = 800;
		if (exactU.empty())
		{
			return;
		}

		Parser exactParser(exactU.c_str());
		Expression exactExpression = exactParser.parse();
		plot.hasExact = true;
		plot.exactX.reserve(exactSamples);
		plot.exactY.reserve(exactSamples);

		for (int i = 0; i < exactSamples; ++i)
		{
			const double x = length * i / (exactSamples - 1);
			plot.exactX.push_back(x);
			plot.exactY.push_back(eval(exactExpression, x));
		}
	}

	PlotSeries buildPlotSeries(
		const ProblemDefinition& problem,
		const std::vector<SolutionPiece>& pieces,
		const std::vector<std::pair<double, double>>& nodes)
	{
		const int samplesPerPiece = 40;
		PlotSeries plot;
		plot.femX.reserve(static_cast<int>(pieces.size()) * samplesPerPiece);
		plot.femY.reserve(static_cast<int>(pieces.size()) * samplesPerPiece);
		plot.nodeX.reserve(nodes.size());
		plot.nodeY.reserve(nodes.size());

		for (const SolutionPiece& piece : pieces)
		{
			for (int local = 0; local < samplesPerPiece; ++local)
			{
				const double t = static_cast<double>(local) / (samplesPerPiece - 1);
				const double x = piece.left + (piece.right - piece.left) * t;
				plot.femX.push_back(x);
				plot.femY.push_back(piece.freeTerm + piece.xCoeff * x);
			}
		}

		for (const auto& node : nodes)
		{
			plot.nodeX.push_back(node.first);
			plot.nodeY.push_back(node.second);
		}

		appendExactSamples(plot, problem.exactU, problem.length);
		return plot;
	}

	double maxAbsPlotValue(const PlotSeries& plot, double limiter)
	{
		double result = std::max(1.0, std::abs(limiter));
		for (double value : plot.femY)
		{
			result = std::max(result, std::abs(value));
		}
		for (double value : plot.exactY)
		{
			result = std::max(result, std::abs(value));
		}
		return result * 1.1;
	}

	ProblemDefinition buildProblem(const GuiState& state)
	{
		ProblemDefinition problem;
		problem.intervals = state.intervals;
		problem.length = state.length;
		problem.limiter = state.limiter;
		problem.p = state.p;
		problem.q = state.q;
		problem.exactU = state.exactU;

		if (state.mode == 2)
		{
			problem.generatedRightSide = true;
			problem.name = "Generated from exact solution";
			problem.f = buildGeneratedRightSideExpression(problem.p, problem.q, problem.exactU);
		}
		else
		{
			problem.name = "Manual input";
			problem.f = state.f;
		}

		return problem;
	}

	void rebuildPlotSeries(GuiState& state)
	{
		std::vector<SolutionPiece> pieces;
		std::vector<std::pair<double, double>> nodes;
		pieces.reserve(state.problem.intervals);
		nodes.reserve(state.problem.intervals + 1);
		for (int interval = 0; interval < state.problem.intervals; ++interval)
		{
			SolutionPiece piece;
			piece.left = state.problem.length * interval / state.problem.intervals;
			piece.right = state.problem.length * (interval + 1) / state.problem.intervals;
			piece.freeTerm = state.result.solutionByInterval[interval][0];
			piece.xCoeff = state.result.solutionByInterval[interval][1];
			pieces.push_back(piece);
		}

		for (int node = 0; node <= state.problem.intervals; ++node)
		{
			const double x = state.problem.length * node / state.problem.intervals;
			const int interval = std::min(node, state.problem.intervals - 1);
			nodes.emplace_back(x, evaluatePiece(state.result.solutionByInterval[interval], x));
		}

		state.plot = buildPlotSeries(state.problem, pieces, nodes);
	}

	LoadedPlot loadPlotFile(const std::string& path)
	{
		std::ifstream input(path);
		if (!input)
		{
			throw std::runtime_error("Cannot open plot file: " + path);
		}

		LoadedPlot loaded;
		loaded.path = path;
		std::vector<SolutionPiece> pieces;
		std::vector<std::pair<double, double>> nodes;
		std::string section = "meta";
		std::string line;

		while (std::getline(input, line))
		{
			line = trim(line);
			if (line.empty() || line[0] == '#')
			{
				continue;
			}
			if (line.front() == '[' && line.back() == ']')
			{
				section = line.substr(1, line.size() - 2);
				continue;
			}

			if (section == "meta")
			{
				const std::size_t separator = line.find(' ');
				if (separator == std::string::npos)
				{
					continue;
				}
				const std::string key = line.substr(0, separator);
				const std::string value = trim(line.substr(separator + 1));

				if (key == "name") loaded.problem.name = value;
				else if (key == "N") loaded.problem.intervals = std::stoi(value);
				else if (key == "l") loaded.problem.length = std::stod(value);
				else if (key == "m") loaded.problem.limiter = std::stod(value);
				else if (key == "p") loaded.problem.p = value;
				else if (key == "q") loaded.problem.q = value;
				else if (key == "f") loaded.problem.f = value;
				else if (key == "exact_u") loaded.problem.exactU = value;
				else if (key == "case") loaded.boundaryCase = value;
			}
			else if (section == "pieces")
			{
				std::istringstream row(line);
				SolutionPiece piece;
				row >> piece.left >> piece.right >> piece.freeTerm >> piece.xCoeff;
				if (!row.fail())
				{
					pieces.push_back(piece);
				}
			}
			else if (section == "nodes")
			{
				std::istringstream row(line);
				double x = 0.0;
				double value = 0.0;
				row >> x >> value;
				if (!row.fail())
				{
					nodes.emplace_back(x, value);
				}
			}
		}

		if (pieces.empty() || nodes.empty())
		{
			throw std::runtime_error("Plot file does not contain graph data");
		}

		loaded.plot = buildPlotSeries(loaded.problem, pieces, nodes);
		loaded.loaded = true;
		loaded.resetPlotView = true;
		return loaded;
	}

	std::string resultSummary(const GuiState& state)
	{
		if (state.result.boundaryCase == "free end")
		{
			return "Solved: no boundary contact.";
		}
		if (state.result.boundaryCase == "upper boundary")
		{
			return "Solved: contact at upper boundary, x = l.";
		}
		if (state.result.boundaryCase == "lower boundary")
		{
			return "Solved: contact at lower boundary, x = l.";
		}
		return "Solved: " + state.result.boundaryCase + ".";
	}

	void saveCurrentPlot(GuiState& state)
	{
		if (!state.hasResult)
		{
			state.status = "Solve the problem before saving a plot.";
			return;
		}

		const std::string path = showSavePlotDialog();
		if (path.empty())
		{
			return;
		}

		writePlotFile(path, state.problem, state.result);
		state.status = "Plot saved: " + path;
	}

	void loadComparisonPlot(GuiState& state, int index)
	{
		const std::string path = showOpenPlotDialog();
		if (path.empty())
		{
			return;
		}

		state.comparisonPlots[index] = loadPlotFile(path);
		state.comparisonStatus = "Loaded plot " + std::to_string(index + 1) + ": " + path;
	}

	void runPythonPlot()
	{
		const int exitCode = std::system("python plot_solution.py");
		if (exitCode != 0)
		{
			std::cout << "Plot script failed. You can run: python plot_solution.py" << std::endl;
		}
	}

	void solveFromGui(GuiState& state)
	{
		state.problem = buildProblem(state);
		state.generatedF = state.problem.generatedRightSide ? state.problem.f : "";
		state.result = solveProblem(state.problem);

		std::ofstream out("output.txt");
		if (!out)
		{
			throw std::runtime_error("Cannot open output.txt");
		}

		writeReport(out, state.problem, state.result);
		writePlotFile("output_for_plot.txt", state.problem, state.result);
		rebuildPlotSeries(state);
		state.hasResult = true;
		state.resetPlotView = true;

		std::cout << resultSummary(state) << std::endl;
		std::cout << "Done. See output.txt, output_for_plot.txt and solution_plot.png" << std::endl;

		if (state.updatePythonPlot)
		{
			runPythonPlot();
		}

		state.status = resultSummary(state);
	}

	void trySolve(GuiState& state)
	{
		try
		{
			solveFromGui(state);
		}
		catch (const std::exception& error)
		{
			state.hasResult = false;
			state.status = std::string("Error: ") + error.what();
			std::cerr << state.status << std::endl;
		}
	}

	bool canBuildGeneratedRightSide(const GuiState& state)
	{
		return !state.p.empty() && !state.q.empty() && !state.exactU.empty();
	}

	void drawPlotContent(
		const char* title,
		const PlotSeries& plot,
		double length,
		double limiter,
		bool& resetPlotView)
	{
		ImPlot::PushStyleColor(ImPlotCol_PlotBg, ImVec4(0.055f, 0.058f, 0.066f, 1.0f));
		ImPlot::PushStyleColor(ImPlotCol_FrameBg, ImVec4(0.055f, 0.058f, 0.066f, 1.0f));
		ImPlot::PushStyleColor(ImPlotCol_AxisGrid, ImVec4(0.42f, 0.44f, 0.48f, 0.34f));
		ImPlot::PushStyleColor(ImPlotCol_AxisText, ImVec4(0.86f, 0.88f, 0.91f, 1.0f));
		ImPlot::PushStyleColor(ImPlotCol_LegendBg, ImVec4(0.085f, 0.090f, 0.100f, 0.88f));
		ImPlot::PushStyleColor(ImPlotCol_LegendBorder, ImVec4(0.35f, 0.37f, 0.41f, 0.70f));

		if (ImPlot::BeginPlot(title, ImVec2(-1, -1), ImPlotFlags_NoMenus))
		{
			ImPlot::SetupAxes("x", "u(x)", ImPlotAxisFlags_None, ImPlotAxisFlags_None);
			const double yLimit = maxAbsPlotValue(plot, limiter);
			if (resetPlotView)
			{
				ImPlot::SetupAxesLimits(0.0, length, -yLimit, yLimit, ImPlotCond_Always);
				resetPlotView = false;
			}

			const double axisX[2] = { 0.0, length };
			const double axisY[2] = { 0.0, 0.0 };
			const double verticalAxisX[2] = { 0.0, 0.0 };
			const double verticalAxisY[2] = { -yLimit, yLimit };
			ImPlot::SetNextLineStyle(ImVec4(0.78f, 0.80f, 0.84f, 0.80f), 1.0f);
			ImPlot::PlotLine("##x-axis", axisX, axisY, 2);
			ImPlot::SetNextLineStyle(ImVec4(0.78f, 0.80f, 0.84f, 0.80f), 1.0f);
			ImPlot::PlotLine("##y-axis", verticalAxisX, verticalAxisY, 2);

			ImPlot::SetNextLineStyle(ImVec4(0.18f, 0.55f, 1.0f, 1.0f), 2.4f);
			ImPlot::PlotLine("u_h(x)", plot.femX.data(), plot.femY.data(), static_cast<int>(plot.femX.size()));

			ImPlot::SetNextMarkerStyle(ImPlotMarker_Circle, 3.5f, ImVec4(0.18f, 0.55f, 1.0f, 1.0f));
			ImPlot::PlotScatter("mesh nodes", plot.nodeX.data(), plot.nodeY.data(), static_cast<int>(plot.nodeX.size()));

			if (plot.hasExact)
			{
				ImPlot::SetNextLineStyle(ImVec4(1.0f, 0.42f, 0.18f, 1.0f), 1.8f);
				ImPlot::PlotLine("exact u(x)", plot.exactX.data(), plot.exactY.data(), static_cast<int>(plot.exactX.size()));
			}

			const double limiterX[2] = { 0.0, length };
			const double upperLimiterY[2] = { limiter, limiter };
			const double lowerLimiterY[2] = { -limiter, -limiter };
			ImPlot::SetNextLineStyle(ImVec4(0.82f, 0.76f, 0.45f, 1.0f), 1.0f);
			ImPlot::PlotLine("limiter", limiterX, upperLimiterY, 2);
			ImPlot::SetNextLineStyle(ImVec4(0.82f, 0.76f, 0.45f, 1.0f), 1.0f);
			ImPlot::PlotLine("##lower_limiter", limiterX, lowerLimiterY, 2);

			ImPlot::EndPlot();
		}

		ImPlot::PopStyleColor(6);
	}

	void drawPlotPanel(GuiState& state, const ImVec2& size)
	{
		ImGui::PushStyleColor(ImGuiCol_ChildBg, ImVec4(0.070f, 0.075f, 0.082f, 1.0f));
		ImGui::BeginChild("PlotPanel", size, false, ImGuiWindowFlags_NoScrollbar);

		if (state.hasResult && !state.plot.femX.empty())
		{
			drawPlotContent("FEM solution##current", state.plot, state.problem.length, state.problem.limiter, state.resetPlotView);
		}
		else
		{
			ImGui::TextWrapped("Press Solve to build the FEM graph.");
		}

		ImGui::EndChild();
		ImGui::PopStyleColor();
	}

	void drawModePanel(GuiState& state, const ImVec2& size)
	{
		ImGui::PushStyleColor(ImGuiCol_ChildBg, ImVec4(0.075f, 0.083f, 0.095f, 1.0f));
		ImGui::BeginChild("ModePanel", size, false);

		ImGui::TextUnformatted("Task");
		ImGui::Separator();
		const char* tasks[] = {
			"1 - p(x), q(x), f(x)",
			"2 - exact u(x), generate f(x)",
			"3 - comparison of graphs"
		};
		int selectedTask = state.mode - 1;
		if (ImGui::Combo("Mode", &selectedTask, tasks, IM_ARRAYSIZE(tasks)))
		{
			state.mode = selectedTask + 1;
		}
		if (state.mode != 3)
		{
			ImGui::Checkbox("Update Python plot after solve", &state.updatePythonPlot);
		}

		ImGui::EndChild();
		ImGui::PopStyleColor();
	}

	void drawInputPanel(GuiState& state, const ImVec2& size)
	{
		ImGui::PushStyleColor(ImGuiCol_ChildBg, ImVec4(0.075f, 0.073f, 0.065f, 1.0f));
		ImGui::PushStyleColor(ImGuiCol_FrameBg, ImVec4(0.125f, 0.122f, 0.105f, 1.0f));
		ImGui::PushStyleColor(ImGuiCol_FrameBgHovered, ImVec4(0.170f, 0.165f, 0.135f, 1.0f));
		ImGui::PushStyleColor(ImGuiCol_FrameBgActive, ImVec4(0.205f, 0.195f, 0.155f, 1.0f));
		ImGui::BeginChild("InputPanel", size, false);

		ImGui::TextUnformatted("Problem data");
		ImGui::Separator();
		ImGui::InputInt("N", &state.intervals);
		ImGui::InputDouble("l", &state.length, 0.0, 0.0, "%.10g");
		ImGui::InputDouble("m", &state.limiter, 0.0, 0.0, "%.10g");
		ImGui::InputText("p(x)", &state.p);
		ImGui::InputText("q(x)", &state.q);

		if (state.mode == 1)
		{
			ImGui::InputText("f(x)", &state.f);
			ImGui::InputText("exact u(x)", &state.exactU);
		}
		else
		{
			ImGui::InputText("exact u(x)", &state.exactU);
			if (canBuildGeneratedRightSide(state))
			{
				try
				{
					const std::string generated = buildGeneratedRightSideDisplayExpression(state.p, state.q, state.exactU);
					ImGui::TextWrapped("f(x) = %s", generated.c_str());
				}
				catch (const std::exception& error)
				{
					ImGui::TextWrapped("f(x) display is unavailable: %s", error.what());
				}
			}
			else
			{
				ImGui::TextWrapped("f(x) will appear after p(x), q(x), and exact u(x) are entered.");
			}
		}

		if (ImGui::Button("Solve", ImVec2(120.0f, 0.0f)))
		{
			trySolve(state);
		}
		ImGui::SameLine();
		if (!state.hasResult)
		{
			ImGui::BeginDisabled();
		}
		if (ImGui::Button("Save plot", ImVec2(120.0f, 0.0f)))
		{
			try
			{
				saveCurrentPlot(state);
			}
			catch (const std::exception& error)
			{
				state.status = std::string("Error: ") + error.what();
			}
		}
		if (!state.hasResult)
		{
			ImGui::EndDisabled();
		}

		ImGui::Separator();
		ImGui::TextWrapped("%s", state.status.c_str());

		ImGui::EndChild();
		ImGui::PopStyleColor(4);
	}

	std::string comparisonSummary(const LoadedPlot& plot)
	{
		if (!plot.loaded)
		{
			return "No plot loaded.";
		}

		std::ostringstream out;
		if (!plot.problem.name.empty())
		{
			out << plot.problem.name << " | ";
		}
		out << "N=" << plot.problem.intervals
			<< ", l=" << plot.problem.length
			<< ", m=" << plot.problem.limiter
			<< " | p=" << plot.problem.p
			<< " | q=" << plot.problem.q;
		if (!plot.problem.f.empty() && plot.problem.f.find("dx(") == std::string::npos)
		{
			out << " | f=" << plot.problem.f;
		}
		if (!plot.problem.exactU.empty())
		{
			out << " | u=" << plot.problem.exactU;
		}
		if (!plot.boundaryCase.empty())
		{
			out << " | " << plot.boundaryCase;
		}
		return out.str();
	}

	void drawLoadedPlotPanel(const char* title, LoadedPlot& loaded, const ImVec2& size)
	{
		ImGui::PushStyleColor(ImGuiCol_ChildBg, ImVec4(0.070f, 0.075f, 0.082f, 1.0f));
		ImGui::BeginChild(title, size, false);

		if (loaded.loaded)
		{
			ImGui::TextWrapped("%s", comparisonSummary(loaded).c_str());
			ImGui::Separator();
			drawPlotContent(
				title,
				loaded.plot,
				loaded.problem.length,
				loaded.problem.limiter,
				loaded.resetPlotView);
		}
		else
		{
			ImGui::TextUnformatted("Load a saved plot.");
		}

		ImGui::EndChild();
		ImGui::PopStyleColor();
	}

	void drawComparisonUi(GuiState& state, const ImVec2& available)
	{
		const float gap = 8.0f;
		const float toolbarHeight = 96.0f;

		ImGui::PushStyleColor(ImGuiCol_ChildBg, ImVec4(0.075f, 0.083f, 0.095f, 1.0f));
		ImGui::BeginChild("ComparisonToolbar", ImVec2(-1.0f, toolbarHeight), false);
		ImGui::TextUnformatted("Comparison of Graphs");
		ImGui::Separator();

		const char* tasks[] = {
			"1 - p(x), q(x), f(x)",
			"2 - exact u(x), generate f(x)",
			"3 - comparison of graphs"
		};
		int selectedTask = state.mode - 1;
		ImGui::SetNextItemWidth(280.0f);
		if (ImGui::Combo("Mode", &selectedTask, tasks, IM_ARRAYSIZE(tasks)))
		{
			state.mode = selectedTask + 1;
		}

		ImGui::SameLine();
		if (ImGui::Button("Load left plot", ImVec2(135.0f, 0.0f)))
		{
			try
			{
				loadComparisonPlot(state, 0);
			}
			catch (const std::exception& error)
			{
				state.comparisonStatus = std::string("Error: ") + error.what();
			}
		}
		ImGui::SameLine();
		if (ImGui::Button("Load right plot", ImVec2(135.0f, 0.0f)))
		{
			try
			{
				loadComparisonPlot(state, 1);
			}
			catch (const std::exception& error)
			{
				state.comparisonStatus = std::string("Error: ") + error.what();
			}
		}

		ImGui::TextWrapped("%s", state.comparisonStatus.c_str());
		ImGui::EndChild();
		ImGui::PopStyleColor();

		const float plotHeight = std::max(160.0f, available.y - toolbarHeight - gap);
		const float splitterWidth = 8.0f;
		const float plotsWidth = std::max(240.0f, available.x - splitterWidth - gap * 2.0f);
		const float leftWidth = std::max(80.0f, plotsWidth * state.comparisonSplit);
		const float rightWidth = std::max(80.0f, plotsWidth - leftWidth);

		drawLoadedPlotPanel("Left graph##comparison_left", state.comparisonPlots[0], ImVec2(leftWidth, plotHeight));
		ImGui::SameLine(0.0f, gap);

		ImGui::InvisibleButton("##comparison_splitter", ImVec2(splitterWidth, plotHeight));
		if (ImGui::IsItemHovered())
		{
			ImGui::SetMouseCursor(ImGuiMouseCursor_ResizeEW);
		}
		if (ImGui::IsItemActive())
		{
			state.comparisonSplit = std::clamp(
				state.comparisonSplit + ImGui::GetIO().MouseDelta.x / plotsWidth,
				0.08f,
				0.92f);
		}

		ImGui::SameLine(0.0f, gap);
		drawLoadedPlotPanel("Right graph##comparison_right", state.comparisonPlots[1], ImVec2(rightWidth, plotHeight));
	}

	void drawMainUi(GuiState& state)
	{
		const ImGuiViewport* viewport = ImGui::GetMainViewport();
		ImGui::SetNextWindowPos(viewport->WorkPos);
		ImGui::SetNextWindowSize(viewport->WorkSize);

		const ImGuiWindowFlags flags =
			ImGuiWindowFlags_NoDecoration |
			ImGuiWindowFlags_NoMove |
			ImGuiWindowFlags_NoSavedSettings |
			ImGuiWindowFlags_NoBringToFrontOnFocus;

		ImGui::PushStyleColor(ImGuiCol_WindowBg, ImVec4(0.0f, 0.0f, 0.0f, 1.0f));
		ImGui::Begin("AoFEMfPoDoSSNC", nullptr, flags);

		const ImVec2 available = ImGui::GetContentRegionAvail();
		if (state.mode == 3)
		{
			drawComparisonUi(state, available);
			ImGui::End();
			ImGui::PopStyleColor();
			return;
		}

		const float gap = 8.0f;
		const float rightWidth = std::max(410.0f, available.x * 0.38f);
		const float leftWidth = std::max(320.0f, available.x - rightWidth - gap);
		const float modeHeight = std::max(135.0f, available.y * 0.18f);

		drawPlotPanel(state, ImVec2(leftWidth, available.y));
		ImGui::SameLine(0.0f, gap);

		ImGui::BeginGroup();
		drawModePanel(state, ImVec2(rightWidth, modeHeight));
		drawInputPanel(state, ImVec2(rightWidth, available.y - modeHeight - gap));
		ImGui::EndGroup();

		ImGui::End();
		ImGui::PopStyleColor();
	}

	LRESULT WINAPI wndProc(HWND hwnd, UINT msg, WPARAM wParam, LPARAM lParam)
	{
		if (ImGui_ImplWin32_WndProcHandler(hwnd, msg, wParam, lParam))
		{
			return true;
		}

		switch (msg)
		{
		case WM_SIZE:
			if (gDevice != nullptr && wParam != SIZE_MINIMIZED)
			{
				cleanupRenderTarget();
				gSwapChain->ResizeBuffers(0, LOWORD(lParam), HIWORD(lParam), DXGI_FORMAT_UNKNOWN, 0);
				createRenderTarget();
			}
			return 0;
		case WM_SYSCOMMAND:
			if ((wParam & 0xfff0) == SC_KEYMENU)
			{
				return 0;
			}
			break;
		case WM_DESTROY:
			PostQuitMessage(0);
			return 0;
		default:
			break;
		}

		return DefWindowProcW(hwnd, msg, wParam, lParam);
	}
}

int runFemGuiApp()
{
	WNDCLASSEXW wc = {};
	wc.cbSize = sizeof(wc);
	wc.style = CS_CLASSDC;
	wc.lpfnWndProc = wndProc;
	wc.hInstance = GetModuleHandleW(nullptr);
	wc.lpszClassName = L"AoFEMfPoDoSSNCWindow";
	RegisterClassExW(&wc);

	HWND hwnd = CreateWindowW(
		wc.lpszClassName,
		L"AoFEMfPoDoSSNC",
		WS_OVERLAPPEDWINDOW,
		100,
		100,
		1600,
		900,
		nullptr,
		nullptr,
		wc.hInstance,
		nullptr);

	if (!createDeviceD3D(hwnd))
	{
		cleanupDeviceD3D();
		UnregisterClassW(wc.lpszClassName, wc.hInstance);
		return 1;
	}

	ShowWindow(hwnd, SW_SHOWDEFAULT);
	UpdateWindow(hwnd);

	IMGUI_CHECKVERSION();
	ImGui::CreateContext();
	ImPlot::CreateContext();
	ImGuiIO& io = ImGui::GetIO();
	io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;

	ImGui::StyleColorsDark();
	ImGuiStyle& style = ImGui::GetStyle();
	style.WindowRounding = 0.0f;
	style.ChildRounding = 0.0f;
	style.FrameRounding = 3.0f;

	ImGui_ImplWin32_Init(hwnd);
	ImGui_ImplDX11_Init(gDevice, gDeviceContext);

	GuiState state;
	trySolve(state);

	bool done = false;
	while (!done)
	{
		MSG msg;
		while (PeekMessageW(&msg, nullptr, 0U, 0U, PM_REMOVE))
		{
			TranslateMessage(&msg);
			DispatchMessageW(&msg);
			if (msg.message == WM_QUIT)
			{
				done = true;
			}
		}
		if (done)
		{
			break;
		}

		ImGui_ImplDX11_NewFrame();
		ImGui_ImplWin32_NewFrame();
		ImGui::NewFrame();

		drawMainUi(state);

		ImGui::Render();
		const float clearColor[4] = { 0.0f, 0.0f, 0.0f, 1.0f };
		gDeviceContext->OMSetRenderTargets(1, &gMainRenderTargetView, nullptr);
		gDeviceContext->ClearRenderTargetView(gMainRenderTargetView, clearColor);
		ImGui_ImplDX11_RenderDrawData(ImGui::GetDrawData());
		gSwapChain->Present(1, 0);
	}

	ImGui_ImplDX11_Shutdown();
	ImGui_ImplWin32_Shutdown();
	ImPlot::DestroyContext();
	ImGui::DestroyContext();

	cleanupDeviceD3D();
	DestroyWindow(hwnd);
	UnregisterClassW(wc.lpszClassName, wc.hInstance);

	return 0;
}
