#pragma once

#include "../Functions/functions.h"
#include "../Operations/operations.h"
#include "General.h"
#include "../Utils/Utils.h"
#include <typeinfo>

using namespace utils;

template <typename F1, typename F2>
class Derivative<functions::Logarithm<F1, F2>> : public functions::Abstract
{
public:
	Derivative(const functions::Logarithm<F1, F2>& f)
		: m_f2(f.m_f2), m_f1(f.m_f1), m_df1(f.m_f1), m_df2(f.m_f2)
	{
	}

	double operator()(double x) override
	{
		double dfx1 = 0;
		double dfx2 = 0;
		double fx1 = 0;
		double fx2 = 0;
		double factor1 = 0;
		double factor2 = 0;

		if constexpr (std::is_pointer<Derivative<F1>>::value)
			dfx1 += (*m_df1)(x);
		else
			dfx1 += m_df1(x);

		if constexpr (std::is_pointer<F1>::value)
			factor1 = log(*(m_f1))(x);
		else
			factor1 = log(m_f1(x));

		if constexpr (std::is_pointer<Derivative<F2>>::value)
			dfx2 += (*m_df2)(x);
		else
			dfx2 += m_df2(x);

		if constexpr (std::is_pointer<F2>::value)
			factor2 += log(*(m_f2))(x);
		else
			factor2 += log(m_f2(x));

		if constexpr (std::is_pointer<F1>::value)
			fx1 += (*m_f1)(x);
		else
			fx1 += m_f1(x);

		if constexpr (std::is_pointer<F2>::value)
			fx2 += (*m_f2)(x);
		else
			fx2 += m_f2(x);

		return ((fx1 * dfx2 * factor1 - dfx1 * fx2 * factor2) / (fx1 * fx2 * factor1 * factor1));
	}


	F1 m_f1;
	F2 m_f2;
	Derivative<F1> m_df1;
	Derivative<F2> m_df2;

	typedef operations::Divide<operations::Subtract<operations::Multiply<operations::Multiply<F1,
		typename Derivative<F2>::Type>, functions::Logarithm<functions::Const, F1>>,
		operations::Multiply<operations::Multiply<typename Derivative<F1>::Type, F2>,
		functions::Logarithm<functions::Const, F2>>>,
		operations::Multiply<operations::Multiply<operations::Multiply<F1, F2>,
		functions::Logarithm<functions::Const, F1>>, functions::Logarithm<functions::Const,	F1>>> Type;
	// typename - C»À¿

	Type expression()
	{
		functions::Const e(2.71828182845904523536);
		return ((m_f1 * m_df2.expression() * Log(e, m_f1) - m_df1.expression() * m_f2 * Log(e, m_f2)) /
			(m_f1 * m_f2 * Log(e, m_f1) * Log(e, m_f1)));
	}

};
