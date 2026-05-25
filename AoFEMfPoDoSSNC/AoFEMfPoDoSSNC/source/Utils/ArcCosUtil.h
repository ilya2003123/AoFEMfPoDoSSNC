#pragma once

#include "../Functions/functions.h"

namespace utils
{
	template<typename F>
	functions::ArcCos<F> ACos(F& f)
	{
		return functions::ArcCos<F>(f);
	}

	double ACos(double x)
	{
		return acos(x);
	}
}
