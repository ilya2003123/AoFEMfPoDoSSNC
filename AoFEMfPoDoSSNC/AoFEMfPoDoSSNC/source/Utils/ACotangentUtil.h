#pragma once

#include "../Functions/functions.h"

namespace utils
{
	template<typename F>
	functions::ACotangent<F> ACtg(F& f)
	{
		return functions::ACotangent<F>(f);
	}

	double ACtg(double x)
	{
		return (1 / atan(x));
	}
}