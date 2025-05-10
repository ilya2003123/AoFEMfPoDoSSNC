#pragma once

#include "../Operations/operations.h"
#include "../Functions/functions.h"
#include "General.h"

// По сути повторение Add.h в папке Derivative, только со своими фишками

template<typename F1, typename F2>
class Derivative<operations::Subtract<F1, F2>> : public functions::Abstract
{
public:
	Derivative(const operations::Subtract<F1, F2>& f)
		: m_df1(f.m_f1), m_df2(f.m_f2)
	{
	}

	double operator()(double x) override
	{
		double dfx1 = 0;
		double dfx2 = 0;

		if constexpr (std::is_pointer<Derivative<F1>>::value)
			dfx1 += (*m_df1)(x);
		else
			dfx1 += m_df1(x);

		if constexpr (std::is_pointer<Derivative<F2>>::value)
			dfx2 += (*m_df2)(x);
		else
			dfx2 += m_df2(x);

		return dfx1 - dfx2;
	}

	Derivative<F1> m_df1;
	Derivative<F2> m_df2;

	typedef operations::Subtract<typename Derivative<F1>::Type, typename Derivative<F2>::Type> Type;
	Type expression()
	{
		return m_df1.expression() - m_df2.expression();
	}

};


template<typename F1>
class Derivative<operations::Subtract<F1, functions::Const>> : public functions::Abstract
{
public:
	Derivative(const operations::Subtract<F1, functions::Const>& f)
		: m_df1(f.m_f1)
	{
	}

	double operator()(double x) override
	{
		double dfx1 = 0;

		if constexpr (std::is_pointer<Derivative<F1>>::value)
			dfx1 += (*m_df1)(x);
		else
			dfx1 += m_df1(x);

		return dfx1;
	}

	Derivative<F1> m_df1;

	typedef typename Derivative<F1>::Type Type;
	Type expression()
	{
		return m_df1.expression();
	}

};

template <typename F2>
class Derivative<operations::Subtract<functions::Const, F2>> : public functions::Abstract
{
public:
	Derivative(const operations::Subtract<functions::Const, F2>& f)
		: m_df2(f.m_f2)
	{
	}
	double operator()(double x) override
	{
		double dfx2 = 0;

		if constexpr (std::is_pointer<Derivative<F2>>::value)
			dfx2 += (*m_df2)(x);
		else
			dfx2 += m_df2(x);

		return (-dfx2);
	}

	Derivative<F2> m_df2;

	typedef typename Derivative<F2>::Type Type;
	Type expression()
	{
		return (functions::Const(-1) * m_df2.expression());
	}
};

template<>
class Derivative<operations::Subtract<functions::Const, functions::Const>> : public functions::Abstract
{
public:
	Derivative(const operations::Subtract<functions::Const, functions::Const>&)
	{
	}

	double operator()(double) override
	{
		return 0;
	}

	typedef typename functions::Const Type;
	Type expression()
	{
		return functions::Const(0);
	}

};

template <>
class Derivative<operations::Subtract<functions::Simple, functions::Simple>> : public functions::Abstract
{
public:
	Derivative(const operations::Subtract<functions::Simple, functions::Simple>&)
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
