#pragma once

#include "../Functions/functions.h"
#include "../Operations/operations.h"
#include "General.h"
#include "../Utils/Utils.h"

// Просто честно берём ручками, считаем производную и пишем сюда, и так в каждом файлике!
// И конечно же мучаемся с типами, потому что это шаблоны :D

template <class F>
class Derivative<functions::Power<F>> : public functions::Abstract
{
public:
	Derivative(const functions::Power<F>& f)
		: m_f(f.m_f), m_n(f.m_n), m_df(f.m_f)
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

		return (m_n * pow(fx, m_n - 1) * dfx);
	}

	F m_f;
	double m_n;
	Derivative<F> m_df;

	typedef operations::Multiply<operations::Multiply<functions::Const, functions::Power<F>>,
		typename Derivative<F>::Type> Type;
	Type expression()
	{
		return ((functions::Const(m_n) * utils::Pow(m_f, m_n - 1)) * m_df.expression());
	}

};

template <>
class Derivative<functions::Power<functions::Const>> : public functions::Abstract
{
public:
	Derivative(const functions::Power<functions::Const>& )
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

template <>
class Derivative<functions::Power<functions::Simple>> : public functions::Abstract
{
public:
	Derivative(const functions::Power<functions::Simple>& f)
		: m_n(f.m_n)
	{
	}

	double operator()(double x) override
	{
		return (m_n * pow(x, m_n - 1));
	}

	double m_n;


	typedef operations::Multiply<functions::Const, functions::Power<functions::Simple>> Type;
	Type expression()
	{
		return operations::Multiply<functions::Const, functions::Power<functions::Simple>>
			(functions::Const(m_n), utils::Pow(utils::X, m_n - 1));
	}

};