#pragma once

#include "../Functions/functions.h"
#include "../Operations/operations.h"
#include "../Utils/Utils.h"
#include "General.h"


template <typename F>
class Derivative<functions::Sinus<F>> : public functions::Abstract
{
public:
	Derivative(const functions::Sinus<F>& f)
		: m_f(f.m_f), m_df(f.m_f)
	{
	}

	double operator()(double x) override
	{
		double dfx = 0;
		double fx = 0;

		if constexpr (std::is_pointer<Derivative<F>>::value)
			dfx += (*m_df)(x);
		else
			dfx += m_df(x);

		if constexpr (std::is_pointer<F>::value)
			fx += (*m_f)(x);
		else
			fx += m_f(x);

		return (cos(fx) * dfx);
	}

	F m_f;
	Derivative<F> m_df;

	typedef operations::Multiply<functions::Cosinus<F>, typename Derivative<F>::Type> Type;
	Type expression()
	{
		return (utils::Cos(m_f) * m_df.expression());
	}

};
