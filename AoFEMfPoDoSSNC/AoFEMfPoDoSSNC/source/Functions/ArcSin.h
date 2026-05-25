#pragma once

#include "Abstract.h"
#include <cmath>

namespace functions
{
	template<typename F>
	class ArcSin : public Abstract
	{
	public:
		using Type = ArcSin<F>;

		ArcSin(const F& f)
			: m_f(f)
		{
		}

		double operator()(double x) override
		{
			return std::asin(evaluate(m_f, x));
		}

		F m_f;
	};
}
