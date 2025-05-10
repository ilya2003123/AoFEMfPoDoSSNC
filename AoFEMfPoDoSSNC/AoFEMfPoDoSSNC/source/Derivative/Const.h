#pragma once

#include"../Functions/functions.h"
#include "General.h"

template<>  // частичная специализация
class Derivative<functions::Const> : public functions::Abstract
{
public:
	typedef functions::Const Type;
	Derivative(functions::Const) {}

	double operator()(double) override
	{
		return 0;
	}

	Type expression()
	{
		return functions::Const(0);
	}
};