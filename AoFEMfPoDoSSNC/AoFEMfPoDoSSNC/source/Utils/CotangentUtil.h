#pragma once

#include "../Functions/functions.h"

namespace utils
{
	template<typename F>
	functions::Cotangent<F> Ctg(F& f)
	{
		return functions::Cotangent<F>(f);
	}

	double Ctg(double x)
	{
		return (1 / tan(x));
	}
}