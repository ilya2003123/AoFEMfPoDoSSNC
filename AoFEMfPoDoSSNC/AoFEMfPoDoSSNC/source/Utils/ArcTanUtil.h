#pragma once

#include "../Functions/functions.h"

namespace utils
{
	template<typename F>
	functions::ArcTan<F> ATg(F& f)
	{
		return functions::ArcTan<F>(f);
	}

	double ATg(double x)
	{
		return atan(x);
	}
}
