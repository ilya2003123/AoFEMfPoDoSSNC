#pragma once

#include <cmath>
#include "Abstract.h"

namespace functions
{
	template<typename F>
	class Tangent : public functions::Abstract
	{
	public:
		typedef Tangent<F> Type;
		Tangent(const F& f)
			:m_f(f)
		{
		}

		double operator()(double x) override
		{
			if constexpr (std::is_pointer<F>::value)
				return tan((*m_f)(x));
			else
				return tan(m_f(x));
		}

		F m_f;

	};
}