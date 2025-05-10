#pragma once

#include "../Functions/functions.h"

namespace utils  // Логарифмы тоже могут быть разные
{
	template<typename F1, typename F2>
	functions::Logarithm<F1, F2> Log(F1& f1, F2& f2)
	{
		return functions::Logarithm<F1, F2>(f1, f2);
	}

	double Log(double base, double x)  // опять же double double double 
	{
		return (log(x) / log(base));
	}
}