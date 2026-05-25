#pragma once

#include "../Functions/functions.h"

namespace utils
{
	template<typename F>
	functions::ArcSin<F> ASin(F& f)
	{
		return functions::ArcSin<F>(f);
	}

	double ASin(double x)
	{
		return asin(x);
	}
}
