#pragma once

#include "../Functions/functions.h"

namespace utils
{
	template<typename F>
	functions::Tangent<F> Tg(F& f)
	{
		return functions::Tangent<F>(f);
	}

	double Tg(double x)
	{
		return tan(x);
	}
}