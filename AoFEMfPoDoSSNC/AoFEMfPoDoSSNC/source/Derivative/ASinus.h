#pragma once

#include "../Operations/operations.h"
#include "../Functions/functions.h"
#include "General.h"
#include "../Utils/Utils.h"

template <typename F>
class Derivative<functions::ASinus<F>> : public functions::Abstract
{
public:
	Derivative(const functions::ASinus<F>& f)
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

		return (dfx / sqrt(1 - fx * fx));
	}

	F m_f;
	Derivative<F> m_df;

	typedef operations::Divide<typename Derivative<F>::Type,
		functions::Power<operations::Subtract<functions::Const, operations::Multiply<F, F>>>> Type;
	Type expression()
	{
		operations::Subtract<functions::Const, operations::Multiply<F, F>> a = (1 - m_f * m_f);
		return (m_df.expression() / utils::Sqrt(a));
	}

};
