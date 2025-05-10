#pragma once

#include <cmath>
#include "Abstract.h"

namespace functions
{
	template<typename F1, typename F2>
	class Exponent_Power : public functions::Abstract
	{
	public:
		typedef Exponent_Power<F1, F2> Type;
		Exponent_Power(const F1& f1, const F2& f2)
			: m_f1(f1), m_f2(f2)
		{
		}

		double operator()(double x) override
		{
			double fx1 = 0;
			double fx2 = 0;

			if constexpr (std::is_pointer<F1>::value)
				fx1 += (*m_f1)(x);
			else
				fx1 += m_f1(x);

			if constexpr (std::is_pointer<F2>::value)
				fx2 += (*m_f2)(x);
			else
				fx2 += m_f2(x);

			return pow(fx1, fx2);
		}


		F1 m_f1;
		F2 m_f2;

	};

}