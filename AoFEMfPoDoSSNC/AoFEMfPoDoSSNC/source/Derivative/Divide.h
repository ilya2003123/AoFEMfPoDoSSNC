#pragma once

#include "../Functions/functions.h"
#include "../Operations/operations.h"
#include "../Utils/Utils.h"
#include "General.h"

template <typename F1, typename F2>
class Derivative<operations::Divide<F1, F2>> : public functions::Abstract
{
public:
	Derivative(const operations::Divide<F1, F2>& f)
		: m_f1(f.m_f1), m_f2(f.m_f2), m_df1(f.m_f1), m_df2(f.m_f2)
	{
	}

	double operator()(double x) override
	{
		double fx1 = 0;
		double fx2 = 0;
		double dfx1 = 0;
		double dfx2 = 0;

		if constexpr (std::is_pointer<F1>::value)
			fx1 += (*m_f1)(x);
		else
			fx1 += m_f1(x);

		if constexpr (std::is_pointer<F2>::value)
			fx2 += (*m_f2)(x);
		else
			fx2 += m_f2(x);

		if constexpr (std::is_pointer<Derivative<F1>>::value)
			dfx1 += (*m_df1)(x);
		else
			dfx1 += m_df1(x);

		if constexpr (std::is_pointer<Derivative<F2>>::value)
			dfx2 += (*m_df2)(x);
		else
			dfx2 += m_df2(x);

		return ((dfx1 * fx2 - fx1 * dfx2) / (fx2 * fx2));
	}

	F1 m_f1;
	F2 m_f2;
	Derivative<F1> m_df1;
	Derivative<F2> m_df2;

	typedef operations::Multiply<operations::Subtract<operations::Multiply<typename Derivative<F1>::Type, F2>,
		operations::Multiply<F1, typename Derivative<F2>::Type>>, functions::Power<F2>> Type;
	Type expression()
	{
		return ((m_df1.expression() * m_f2 - m_f1 * m_df2.expression()) * utils::Pow(m_f2, -2));
	}
};

template <typename F2>
class Derivative<operations::Divide<functions::Const, F2>> : public functions::Abstract
{
public:
	Derivative(const operations::Divide<functions::Const, F2>& f)
		: m_f1(f.m_f1), m_f2(f.m_f2), m_df2(f.m_f2)
	{
	}

	double operator()(double x) override
	{
		double cfx1 = 0;
		double fx2 = 0;
		double dfx2 = 0;

		if constexpr (std::is_pointer<functions::Const>::value)// По сути не нужно :D Const всегда НЕ указатель
			cfx1 += (*m_f1.m_const);
		else
			cfx1 += m_f1.m_const;

		if constexpr (std::is_pointer<F2>::value)
			fx2 += (*m_f2)(x);
		else
			fx2 += m_f2(x);

		if constexpr (std::is_pointer<Derivative<F2>>::value)
			dfx2 += (*m_df2)(x);
		else
			dfx2 += m_df2(x);

		return ((-cfx1 * dfx2) / (fx2 * fx2));
	}

	functions::Const m_f1;
	F2 m_f2;
	Derivative<F2> m_df2;

	typedef operations::Multiply<functions::Const, operations::Multiply<typename Derivative<F2>::Type,
		functions::Power<F2>>> Type;
	Type expression()
	{
		return (functions::Const(-m_f1.m_const) * (m_df2.expression() * utils::Pow(m_f2, -2)));
	}
};

template <typename F1>
class Derivative<operations::Divide<F1, functions::Const>> : public functions::Abstract
{
public:
	Derivative(const operations::Divide<F1, functions::Const>& f)
		: m_f2(1 / f.m_f2.m_const), m_df1(f.m_f1)
	{
	}

	double operator()(double x) override
	{
		double cfx2 = 0;
		double dfx1 = 0;

		if constexpr (std::is_pointer<functions::Const>::value)
			cfx2 += (*m_f2.m_const);
		else
			cfx2 += m_f2.m_const;

		if constexpr (std::is_pointer<Derivative<F1>>::value)
			dfx1 += (*m_df1)(x);
		else
			dfx1 += m_df1(x);

		return cfx2 * dfx1;
	}

	functions::Const m_f2;
	Derivative<F1> m_df1;

	typedef operations::Multiply<functions::Const, typename Derivative<F1>::Type> Type;
	Type expression()
	{
		return (m_f2 * m_df1.expression());
	}
};

template <>
class Derivative<operations::Divide<functions::Const, functions::Const>> : public functions::Abstract
{
public:
	Derivative(const operations::Divide<functions::Const, functions::Const>& f)
		: m_const(f.m_f1.m_const / f.m_f2.m_const)
	{
	}

	double operator()(double) override
	{
		return 0;
	}

	functions::Const m_const;

	typedef functions::Const Type;
	Type expression()
	{
		return functions::Const(0);
	}
};