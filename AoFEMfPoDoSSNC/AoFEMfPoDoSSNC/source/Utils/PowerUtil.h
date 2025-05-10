#pragma once

#include "../Functions/functions.h"

namespace utils
{
	template<typename F, typename Num>
	functions::Power<F> Pow(F& f, Num n)
	{
		return functions::Power<F>(f, n);
	}

	double Pow(double x, double y)
	{
		return pow(x, y);
	}
}