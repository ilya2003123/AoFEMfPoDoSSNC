#pragma once

#include "../Operations/operations.h"

namespace utils
{
	template<typename F> // Шаблонная функция возведения в квадрат
	operations::Multiply<F, F> Sqr(F& f)
	{
		return operations::Multiply<F, F>(f, f);
	}

	double Sqr(double x)  // double, потому что double чаще всего будет передаваться, перегрузил просто для скорости
	{
		return x * x;
	}

}