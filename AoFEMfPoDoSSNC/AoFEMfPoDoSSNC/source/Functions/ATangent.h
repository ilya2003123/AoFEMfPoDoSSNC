#pragma once

#include <cmath>
#include "Abstract.h"

namespace functions
{
	template<typename F>
	class ATangent : public functions::Abstract
	{
	public:
		typedef ATangent<F> Type;
		ATangent(const F& f)
			:m_f(f)
		{
		}

		double operator()(double x) override
		{
			if constexpr (std::is_pointer<F>::value)
				return atan((*m_f)(x));
			else
				return atan(m_f(x));
		}

		F m_f;

	};
}