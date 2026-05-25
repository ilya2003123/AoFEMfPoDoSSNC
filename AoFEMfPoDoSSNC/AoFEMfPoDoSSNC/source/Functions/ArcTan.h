#pragma once

#include "Abstract.h"
#include <cmath>

namespace functions
{
	template<typename F>
	class ArcTan : public Abstract
	{
	public:
		using Type = ArcTan<F>;

		ArcTan(const F& f)
			: m_f(f)
		{
		}

		double operator()(double x) override
		{
			return std::atan(evaluate(m_f, x));
		}

		F m_f;
	};
}
