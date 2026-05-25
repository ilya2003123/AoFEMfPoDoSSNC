#pragma once

#include "Abstract.h"
#include <cmath>

namespace functions
{
	template<typename F>
	class Tan : public Abstract
	{
	public:
		using Type = Tan<F>;

		Tan(const F& f)
			: m_f(f)
		{
		}

		double operator()(double x) override
		{
			return std::tan(evaluate(m_f, x));
		}

		F m_f;
	};
}
