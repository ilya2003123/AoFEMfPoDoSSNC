#pragma once

#include "Abstract.h"
#include <cmath>

namespace functions
{
	template<typename F1, typename F2>
	class Logarithm : public Abstract
	{
	public:
		using Type = Logarithm<F1, F2>;

		Logarithm(const F1& base, const F2& argument)
			: m_f1(base), m_f2(argument)
		{
		}

		double operator()(double x) override
		{
			return std::log(evaluate(m_f2, x)) / std::log(evaluate(m_f1, x));
		}

		F1 m_f1;
		F2 m_f2;
	};
}
