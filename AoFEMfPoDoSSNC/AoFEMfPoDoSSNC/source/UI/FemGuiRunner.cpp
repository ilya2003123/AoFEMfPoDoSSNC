#define NOMINMAX

#include "FemGuiApp.h"

#include "FemGuiState.h"
#include "imgui.h"
#include "imgui_impl_dx11.h"
#include "imgui_impl_win32.h"
#include "implot.h"

#include <d3d11.h>
#include <windows.h>

extern IMGUI_IMPL_API LRESULT ImGui_ImplWin32_WndProcHandler(HWND hwnd, UINT msg, WPARAM wParam, LPARAM lParam);

namespace
{
	ID3D11Device* gDevice = nullptr;
	ID3D11DeviceContext* gDeviceContext = nullptr;
	IDXGISwapChain* gSwapChain = nullptr;
	ID3D11RenderTargetView* gMainRenderTargetView = nullptr;

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
	initializeGuiState(state);

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

		drawFemGui(state);

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
