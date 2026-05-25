#pragma once

#include "Abstract.h"
#include <cmath>

namespace functions
{
	template<typename F>
	class ArcCot : public Abstract
	{
	public:
		using Type = ArcCot<F>;

		ArcCot(const F& f)
			: m_f(f)
		{
		}

		double operator()(double x) override
		{
			return std::atan(1.0 / evaluate(m_f, x));
		}

		F m_f;
	};
}
