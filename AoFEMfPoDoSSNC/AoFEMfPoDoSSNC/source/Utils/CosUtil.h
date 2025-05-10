#pragma once

#include "../Functions/functions.h"

namespace utils
{
	template<typename F>
	functions::Cosinus<F> Cos(F& f)
	{
		return functions::Cosinus<F>(f);
	}

	double Cos(double x)
	{
		return cos(x);
	}
}