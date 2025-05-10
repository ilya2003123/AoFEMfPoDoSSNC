#pragma once

#include "../Operations/operations.h"
#include "../Functions/functions.h"
#include "General.h"
#include "../Utils/Utils.h"

template <typename F>
class Derivative<functions::Cotangent<F>> : public functions::Abstract
{
public:
	Derivative(const functions::Cotangent<F>& f)
		: m_f(f.m_f), m_df(f.m_f)
	{
	}

	double operator()(double x) override
	{
		double fx = 0;
		double dfx = 0;

		if constexpr (std::is_pointer<F>::value)
			fx += (*m_f)(x);
		else
			fx += m_f(x);

		if constexpr (std::is_pointer<Derivative<F>>::value)
			dfx += (*m_df)(x);
		else
			dfx += m_df(x);

		return (-dfx / (sin(fx) * sin(fx)));
	}

	F m_f;
	Derivative<F> m_df;

	typedef operations::Divide<operations::Multiply<functions::Const, typename Derivative<F>::Type>,
		operations::Multiply<functions::Sinus<F>, functions::Sinus<F>>> Type;
	Type expression()
	{
		return (functions::Const(-1) * m_df.expression() / (utils::Sin(m_f) * utils::Sin(m_f)));
	}

};
