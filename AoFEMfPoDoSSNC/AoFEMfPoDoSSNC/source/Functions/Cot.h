#pragma once

#include "Abstract.h"
#include <cmath>

namespace functions
{
	template<typename F>
	class Cot : public Abstract
	{
	public:
		using Type = Cot<F>;

		Cot(const F& f)
			: m_f(f)
		{
		}

		double operator()(double x) override
		{
			return 1.0 / std::tan(evaluate(m_f, x));
		}

		F m_f;
	};
}
