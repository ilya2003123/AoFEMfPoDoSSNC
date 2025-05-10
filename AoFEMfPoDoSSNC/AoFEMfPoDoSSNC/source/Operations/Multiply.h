#pragma once

#include "../Functions/Abstract.h"

namespace operations  // тут тем более
{
	template<typename F1, typename F2>
	class Multiply : public functions::Abstract
	{
	public:
		typedef Multiply<F1, F2> Type;
		Multiply(const F1& f1, const F2& f2)
			:m_f1(f1), m_f2(f2)
		{
		}

		double operator()(double x) override
		{
			double f1 = 0;
			double f2 = 0;

			if constexpr (std::is_pointer<F1>::value)
				f1 += (*m_f1)(x);
			else
				f1 += m_f1(x);

			if constexpr (std::is_pointer_v<F2>)
				f2 += (*m_f2)(x);
			else
				f2 += m_f2(x);

			return f1 * f2;
		}

		F1 m_f1;
		F2 m_f2;

	};
}