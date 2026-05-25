#pragma warning(disable: 26451)

#include "UI/FemGuiApp.h"
#include <clocale>

int main()
{
	setlocale(LC_ALL, "rus");
	setlocale(LC_NUMERIC, "C");
	return runFemGuiApp();
}
