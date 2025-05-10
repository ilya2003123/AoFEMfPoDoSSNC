#pragma once

#include "../Functions/functions.h"
#include "../Operations/operations.h"
#include "General.h"

template<typename F1, typename F2>
class Derivative<operations::Multiply<F1, F2>> : public functions::Abstract
{
public:
	Derivative(const operations::Multiply<F1, F2>& f)
		: m_df1(f.m_f1), m_df2(f.m_f2), m_f1(f.m_f1), m_f2(f.m_f2)
	{
	}

	double operator()(double x) override
	{
		double fx1 = 0;
		double dfx1 = 0;
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

		if constexpr (std::is_pointer<Derivative<F1>>::value)
			dfx1 += (*m_df1)(x);
		else
			dfx1 += m_df1(x);

		if constexpr (std::is_pointer<Derivative<F2>>::value)
			dfx2 += (*m_df2)(x);
		else
			dfx2 += m_df2(x);

		return ((dfx1 * fx2) + (fx1 * dfx2));
	}

	Derivative<F1> m_df1;
	Derivative<F2> m_df2;
	F1 m_f1;
	F2 m_f2;

	typedef operations::Add<operations::Multiply<typename Derivative<F1>::Type, F2>,
		operations::Multiply<F1, typename Derivative <F2>::Type>> Type;
	Type expression()
	{
		return ((m_df1.expression() * m_f2) + (m_f1 * m_df2.expression()));
	}

};


template<typename F1>
class Derivative<operations::Multiply<F1, functions::Const>> : public functions::Abstract
{
public:
	Derivative(const operations::Multiply<F1, functions::Const>& f)
		: m_df1(f.m_f1), m_f2(f.m_f2)
	{
	}

	double operator()(double x) override
	{
		double dfx1 = 0;
		double cfx2 = 0;

		if constexpr (std::is_pointer<functions::Const>::value)
			cfx2 += (*m_f2.m_const);
		else
			cfx2 += m_f2.m_const;

		if constexpr (std::is_pointer<Derivative<F1>>::value)
			dfx1 += (*m_df1)(x);
		else
			dfx1 += m_df1(x);

		return (dfx1 * cfx2);
	}

	Derivative<F1> m_df1;
	functions::Const m_f2;

	typedef operations::Multiply<functions::Const, typename Derivative<F1>::Type> Type;
	Type expression()
	{
		return (m_f2 * m_df1.expression());
	}

};

template<typename F2>
class Derivative<operations::Multiply<functions::Const, F2>> : public functions::Abstract
{
public:

	Derivative(const operations::Multiply<functions::Const, F2>& f)
		: m_f1(f.m_f1), m_df2(f.m_f2)
	{
	}

	double operator()(double x) override
	{
		double cfx1 = 0;
		double dfx2 = 0;

		if constexpr (std::is_pointer<functions::Const>::value)
			cfx1 += (*m_f1.m_const);
		else
			cfx1 += m_f1.m_const;

		if constexpr (std::is_pointer<Derivative<F2>>::value)
			dfx2 += (*m_df2)(x);
		else
			dfx2 += m_df2(x);

		return (cfx1 * dfx2);
	}


	functions::Const m_f1;
	Derivative<F2> m_df2;

	typedef operations::Multiply<functions::Const, typename Derivative<F2>::Type> Type;
	Type expression()
	{
		return (m_f1 * m_df2.expression());
	}

};

template<>
class Derivative<operations::Multiply<functions::Const, functions::Const>> : public functions::Abstract
{
public:
	Derivative(const operations::Multiply<functions::Const, functions::Const>&)
	{
	}

	double operator()(double) override
	{
		return 0;
	}

	typedef functions::Const Type;
	Type expression()
	{
		return functions::Const(0);
	}

};