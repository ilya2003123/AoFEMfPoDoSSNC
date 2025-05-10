#pragma once

#include <cmath>
#include "Abstract.h"

namespace functions
{
	template<typename F>
	class ACosinus : public functions::Abstract
	{
	public:
		typedef ACosinus<F> Type;
		ACosinus(const F& f)
			:m_f(f)
		{
		}

		double operator()(double x) override
		{
			if constexpr (std::is_pointer<F>::value)
				return acos((*m_f)(x));
			else
				return acos(m_f(x));
		}

		F m_f;

	};
}