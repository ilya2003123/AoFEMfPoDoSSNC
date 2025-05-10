#pragma once

#include "../Operations/operations.h"
#include "../Functions/functions.h"
#include "General.h"
#include "../Utils/Utils.h"

template <typename F>
class Derivative<functions::Tangent<F>> : public functions::Abstract
{
public:
	Derivative(const functions::Tangent<F>& f)
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

		return (dfx / (cos(fx) * cos(fx)));
	}

	F m_f;
	Derivative<F> m_df;

	typedef operations::Divide<typename Derivative<F>::Type, operations::Multiply<functions::Cosinus<F>,
		functions::Cosinus<F>>> Type;
	Type expression()
	{
		return (m_df.expression() / (utils::Cos(m_f) * utils::Cos(m_f)));
	}

};