#pragma once

#include "../Functions/functions.h"

namespace utils
{
	template<typename F>
	functions::Power<F> Sqrt(F& f)  // тоже шаблонная функция, извлечения корня
	{
		return functions::Power<F>(f, 0.5);
	}

	double Sqrt(double x)  // Опять же для скорости выполнения программы
	{
		return sqrt(x);
	}

}