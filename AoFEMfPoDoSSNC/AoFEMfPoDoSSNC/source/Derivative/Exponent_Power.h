#pragma once

#include "../Functions/functions.h"
#include "../Operations/operations.h"
#include "General.h"
#include "../Utils/Utils.h"

template <typename F1, typename F2>
class Derivative<functions::Exponent_Power<F1, F2>> : public functions::Abstract
{
public:
	Derivative(const functions::Exponent_Power<F1, F2>& f)
		: m_f1(f.m_f1), m_factor(log(f.m_f1))
		, m_f2(f.m_f2), m_df2(f.m_f2)
	{
	}

	double operator()(double x) override
	{
		double fx1 = 0;
		double fx2 = 0;
		double dfx2 = 0;

		if constexpr (std::is_pointer<F1>::value)
			fx1 += (*m_f1)(x);
		else
			fx1 += m_f1(x);

		if constexpr (std::is_pointer<F2>::value)
			fx2 += (*m_f2)(x);
		else
			fx2 += m_f2(x);

		if constexpr (std::is_pointer<Derivative<F2>>::value)
			dfx2 += (*m_df2)(x);
		else
			dfx2 += m_df2(x);

		return (pow(fx1, fx2) * (m_factor * dfx2 + fx2 / fx1));
	}



	F1 m_f1;
	F2 m_f2;
	Derivative<F2> m_df2;
	double m_factor;

	typedef operations::Multiply<functions::Exponent_Power<F1, F2>,
		operations::Add<operations::Multiply<typename Derivative<F2>::Type, functions::Const>,
		operations::Divide<F2, F1>>> Type;
	Type expression()
	{
		return (utils::Exp_Pow(m_f1, m_f2) * (m_df2.expression() * m_factor + m_f2 / m_f1));
	}


};