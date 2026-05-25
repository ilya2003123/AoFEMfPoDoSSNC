#pragma once

#include "../Functions/Abstract.h"

namespace operations
{
	template<typename F1, typename F2>
	class Divide : public functions::Abstract
	{
	public:
		using Type = Divide<F1, F2>;

		Divide(const F1& f1, const F2& f2)
			: m_f1(f1), m_f2(f2)
		{
		}

		double operator()(double x) override
		{
			return functions::evaluate(m_f1, x) / functions::evaluate(m_f2, x);
		}

		F1 m_f1;
		F2 m_f2;
	};
}
