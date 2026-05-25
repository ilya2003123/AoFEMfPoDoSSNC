#pragma once

#include "Abstract.h"
#include <cmath>

namespace functions
{
	template<typename F>
	class ArcCos : public Abstract
	{
	public:
		using Type = ArcCos<F>;

		ArcCos(const F& f)
			: m_f(f)
		{
		}

		double operator()(double x) override
		{
			return std::acos(evaluate(m_f, x));
		}

		F m_f;
	};
}
