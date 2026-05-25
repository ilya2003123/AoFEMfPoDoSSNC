#pragma once

#include "Abstract.h"
#include <cmath>

namespace functions
{
	template<typename F>
	class Exponent : public Abstract
	{
	public:
		using Type = Exponent<F>;

		Exponent(double base, const F& f)
			: m_base(base), m_f(f)
		{
		}

		double operator()(double x) override
		{
			return std::pow(m_base, evaluate(m_f, x));
		}

	private:
		double m_base;
		F m_f;
	};
}
