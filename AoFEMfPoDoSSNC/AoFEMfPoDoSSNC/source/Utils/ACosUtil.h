#pragma once

#include "../Functions/functions.h"

namespace utils
{
	template<typename F>
	functions::ACosinus<F> ACos(F& f)
	{
		return functions::ACosinus<F>(f);
	}

	double ACos(double x)
	{
		return acos(x);
	}
}