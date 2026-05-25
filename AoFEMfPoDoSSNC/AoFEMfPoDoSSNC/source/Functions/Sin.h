#pragma once

#include "Abstract.h"
#include <cmath>

namespace functions
{
	template<typename F>
	class Sin : public Abstract
	{
	public:
		using Type = Sin<F>;

		Sin(const F& f)
			: m_f(f)
		{
		}

		double operator()(double x) override
		{
			return std::sin(evaluate(m_f, x));
		}

		F m_f;
	};
}
