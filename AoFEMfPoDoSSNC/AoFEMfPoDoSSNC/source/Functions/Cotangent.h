#pragma once

#include <cmath>
#include "Abstract.h"

namespace functions
{
	template<typename F>
	class Cotangent : public functions::Abstract
	{
	public:
		typedef Cotangent<F> Type;
		Cotangent(const F& f)
			:m_f(f)
		{
		}

		double operator()(double x) override
		{
			if constexpr (std::is_pointer<F>::value)
				return (1 / tan((*m_f)(x)));
			else
				return (1 / tan(m_f(x)));
		}

		F m_f;

	};
}