#pragma once

#include "../Functions/functions.h"

namespace utils
{
	template<typename F>
	functions::ATangent<F> ATg(F& f)
	{
		return functions::ATangent<F>(f);
	}

	double ATg(double x)
	{
		return atan(x);
	}
}