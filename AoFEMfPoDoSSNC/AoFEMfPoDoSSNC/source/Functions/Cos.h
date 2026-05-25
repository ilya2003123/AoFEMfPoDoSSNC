#pragma once

#include "Abstract.h"
#include <cmath>

namespace functions
{
	template<typename F>
	class Cos : public Abstract
	{
	public:
		using Type = Cos<F>;

		Cos(const F& f)
			: m_f(f)
		{
		}

		double operator()(double x) override
		{
			return std::cos(evaluate(m_f, x));
		}

		F m_f;
	};
}
