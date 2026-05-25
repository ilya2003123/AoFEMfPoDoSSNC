#pragma once

#include "../Functions/Abstract.h"

namespace operations
{
	template<typename F>
	class UnarMinus : public functions::Abstract
	{
	public:
		using Type = UnarMinus<F>;

		UnarMinus(const F& f)
			: m_f(f)
		{
		}

		double operator()(double x) override
		{
			return -functions::evaluate(m_f, x);
		}

		F m_f;
	};
}
