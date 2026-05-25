#pragma once

#include "../Functions/functions.h"

namespace utils
{
	template<typename F>
	functions::Cos<F> Cos(F& f)
	{
		return functions::Cos<F>(f);
	}

	double Cos(double x)
	{
		return cos(x);
	}
}
