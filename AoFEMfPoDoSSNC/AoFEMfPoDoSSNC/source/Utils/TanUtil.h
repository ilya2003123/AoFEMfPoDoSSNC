#pragma once

#include "../Functions/functions.h"

namespace utils
{
	template<typename F>
	functions::Tan<F> Tg(F& f)
	{
		return functions::Tan<F>(f);
	}

	double Tg(double x)
	{
		return tan(x);
	}
}
