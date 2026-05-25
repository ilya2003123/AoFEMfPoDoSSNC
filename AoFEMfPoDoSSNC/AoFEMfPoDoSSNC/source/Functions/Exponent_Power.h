#pragma once

#include "Abstract.h"
#include <cmath>

namespace functions
{
	template<typename F1, typename F2>
	class Exponent_Power : public Abstract
	{
	public:
		using Type = Exponent_Power<F1, F2>;

		Exponent_Power(const F1& f1, const F2& f2)
			: m_f1(f1), m_f2(f2)
		{
		}

		double operator()(double x) override
		{
			return std::pow(evaluate(m_f1, x), evaluate(m_f2, x));
		}

	private:
		F1 m_f1;
		F2 m_f2;
	};
}
