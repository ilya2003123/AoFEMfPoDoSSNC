#pragma once

#include "../Functions/functions.h"

namespace utils  // показательная функция бывает разная
{
	template<typename F>
	functions::Exponent<F> Exp(double base, F& f)  // Шаблонная функция для double
	{
		return functions::Exponent<F>(base, f);
	}

	template<typename F>
	functions::Exponent<F> Exp(int base, F& f)  // шаблонная функция для int 
	{
		return functions::Exponent<F>(base, f);
	}

	double Exp(double base, double x)  // Если кратко, три слова: double double double
	{
		return pow(base, x);
	}


}