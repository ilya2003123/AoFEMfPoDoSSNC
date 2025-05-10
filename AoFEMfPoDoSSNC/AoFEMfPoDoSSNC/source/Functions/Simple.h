#pragma once

#include "Abstract.h"

namespace functions  // Тот же самый namespace
{
	class Simple : public functions::Abstract
	{
	public:
		typedef Simple Type;   // Тоже обзываю, будет тёской Const'a

		double operator()(double x) override  // Так же перегружаю оператор вызова функции 
		{                                  // double передаю и double получаю, обязательно всё const!!!
			return x;
		}
	};
}