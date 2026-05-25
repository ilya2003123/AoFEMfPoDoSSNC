#pragma once

#include "../Functions/functions.h"

namespace utils
{
	template<typename F>
	functions::ArcCot<F> ACtg(F& f)
	{
		return functions::ArcCot<F>(f);
	}

	double ACtg(double x)
	{
		return (1 / atan(x));
	}
}
