#pragma once

#include "../Functions/functions.h"

namespace utils
{
	template<typename F>
	functions::Sinus<F> Sin(F& f)
	{
		return functions::Sinus<F>(f);
	}

	double Sin(double x)
	{
		return sin(x);
	}
}