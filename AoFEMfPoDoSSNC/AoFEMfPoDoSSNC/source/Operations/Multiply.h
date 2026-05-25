#pragma once

#include "../Functions/Abstract.h"

namespace operations
{
	template<typename F1, typename F2>
	class Multiply : public functions::Abstract
	{
	public:
		using Type = Multiply<F1, F2>;

		Multiply(const F1& f1, const F2& f2)
			: m_f1(f1), m_f2(f2)
		{
		}

		double operator()(double x) override
		{
			return functions::evaluate(m_f1, x) * functions::evaluate(m_f2, x);
		}

		F1 m_f1;
		F2 m_f2;
	};
}
