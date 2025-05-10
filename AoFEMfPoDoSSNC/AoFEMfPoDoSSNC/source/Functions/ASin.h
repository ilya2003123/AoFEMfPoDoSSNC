#pragma once

#include <cmath>
#include "Abstract.h"

namespace functions
{
	template<typename F>
	class ASinus : public functions::Abstract
	{
	public:
		typedef ASinus<F> Type;
		ASinus(const F& f)
			:m_f(f)
		{
		}

		double operator()(double x) override
		{
			if constexpr (std::is_pointer<F>::value)
				return asin((*m_f)(x));
			else
				return asin(m_f(x));
		}

		F m_f;

	};
}