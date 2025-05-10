#pragma once

#include "../Functions/functions.h"
#include "../Operations/operations.h"
#include "General.h"
#include "../Utils/Utils.h"

template <typename F>
class Derivative<functions::Exponent<F>> : public functions::Abstract
{
public:
	Derivative(const functions::Exponent<F>& f)
		: m_base(f.m_base), m_factor(log(f.m_base))
		, m_f(f.m_f), m_df(f.m_f)
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

		return (pow(m_base, fx) * m_factor * dfx);
	}

	double m_base, m_factor;
	F m_f;
	Derivative<F> m_df;

	typedef operations::Multiply<operations::Multiply<functions::Exponent<F>, functions::Const>,
		typename Derivative<F>::Type> Type;
	Type expression() 
	{
		return (utils::Exp(m_base, m_f) * m_factor * m_df.expression());
	}

};