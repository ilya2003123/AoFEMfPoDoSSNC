#pragma once

#include "../Functions/functions.h"

namespace utils
{
	template<typename F1, typename F2>
	functions::Exponent_Power<F1, F2> Exp_Pow(F1& f1, F2& f2)
	{
		return functions::Exponent_Power<F1, F2>(f1, f2);
	}

	double Exp_Pow(double base, double x) 
	{
		return pow(base, x);
	}

}