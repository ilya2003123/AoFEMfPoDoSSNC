#pragma once

#include "../Functions/functions.h"

namespace utils
{
	template<typename F>
	functions::Cot<F> Ctg(F& f)
	{
		return functions::Cot<F>(f);
	}

	double Ctg(double x)
	{
		return (1 / tan(x));
	}
}
