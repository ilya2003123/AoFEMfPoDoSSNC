#pragma once

#include "../Functions/functions.h"

namespace utils
{
	template<typename F>
	functions::Sin<F> Sin(F& f)
	{
		return functions::Sin<F>(f);
	}

	double Sin(double x)
	{
		return sin(x);
	}
}
