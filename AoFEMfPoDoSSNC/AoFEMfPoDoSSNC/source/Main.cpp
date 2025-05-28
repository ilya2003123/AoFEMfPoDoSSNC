#pragma warning(disable: 26451)

#include <iostream>
#include <string>
#include <vector>
#include "Functions/functions.h"
#include "Operations/operations.h"
#include "Utils/Utils.h"
#include "Derivative/derivative.h"
#include "Parser/Parser.h"
#include "test/test.h"
#include "basis/basis.h"
#include "coefficient/coefficient.h"
//#include <ionanip> - для вывода, можно отрегулировать сколько знаков после запятой выводится!

using namespace utils;


int main()
{
	//test();
	setlocale(LC_ALL, "rus");

	std::cout << "Введите колличество разбиений: ";
	int n;
	std::cin >> n;

	std::cout << "Введите границу m: ";
	double m;
	std::cin >> m;

	Cap** phi = new Cap * [n];
	for (int i = 0; i < n; i++)
	{
		phi[i] = new Cap[n + 1];
	}

	equationStraightLine(phi, n);

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n + 1; j++)
		{
			std::cout << phi[i][j];
		}
		std::endl(std::cout);
	}

	// сейчас работает для конретного уравнения кусочной
   // функции, сделай так, чтобы она работала по полной, потом добавь возможность интервала, после чего сделай ввод функций p и q
   // после этого посмотри в гпт, и сделай метод прогонки, потом реализуй всё тоже самое для d_i, и у тебя получится нужное уравнение.

   // или просто глянь в гпт, он там чот умное написал. :D



	std::vector<double> coeffC = thomasAlgorithm(integrateProduct(phi, 1, 1, n), rightSystemCoefficientC(n, 1, phi));
	std::vector<double> coeffD = thomasAlgorithm(integrateProduct(phi, 1, 1, n), rightSystemCoefficientD(n));

	double c = 0.0;

	if (coeffC[n - 1] >= m)
	{
		c = (m - coeffC[n - 1]) / coeffD[n - 1];
	}
	else if (coeffC[n - 1] <= -m)
	{
		c = (-m - coeffC[n - 1]) / coeffD[n - 1];
	}
	else
	{
		c = 0.0;
	}

	std::vector<std::vector<std::vector<double>>> u_h;
	u_h.resize(n, std::vector<std::vector<double>>(n, std::vector<double>(2, 0.0)));

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			u_h[i][j][0] += (coeffC[i] + c * coeffD[i]) * phi[i][j].m_equation.m_coeff[0];
			u_h[i][j][1] += (coeffC[i] + c * coeffD[i]) * phi[i][j].m_equation.m_coeff[1];
		}
	}

	for (int i = 0; i < n; i++)
	{
		std::cout << "System №" << i + 1 << ":\n";
		for (int j = 0; j < n; j++)
		{
			std::cout << "Equation №" << j + 1 << ": " << u_h[i][j][1] << "+" << u_h[i][j][0] << "\n";
		}
	}

	return 0;
}