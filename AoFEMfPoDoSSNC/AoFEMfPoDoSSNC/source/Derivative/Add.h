#pragma once

#include "../Operations/operations.h"
#include "../Functions/functions.h"
#include "General.h"

template<typename F1, typename F2> // Передаём сюда два типа
class Derivative<operations::Add<F1, F2>> : public functions::Abstract // класс от операции сложения двух типов
{
public:
	Derivative(const operations::Add<F1, F2>& f) // конструктор, тут не минималистичен
		: m_df1(f.m_f1), m_df2(f.m_f2)
	{
	}

	double operator()(double x) override
	{
		double df1 = 0;
		double df2 = 0;

		if constexpr (std::is_pointer<Derivative<F1>>::value)
			df1 += (*m_df1)(x);
		else
			df1 += m_df1(x);

		if constexpr (std::is_pointer<Derivative<F2>>::value)
			df2 += (*m_df2)(x);
		else
			df2 += m_df2(x);

		return df1 + df2;
	}

	Derivative<F1> m_df1;
	Derivative<F2> m_df2;

	// Я думаю пришло время рассказать ещё одну вещь, для чего нам везде писать слово typename?
	// Что за слово такое, а вот оно нам надо, чтобы компилятор понимал, что я использую шаблоны
	// т.к. при "прочтении" он "<", ">" представляет как знак меньше или больше соответственно

	typedef operations::Add<typename Derivative<F1>::Type, typename Derivative<F2>::Type> Type;
	Type expression()
	{
		return m_df1.expression() + m_df2.expression();
	}
	// производная от суммы "функций" будет сумма производных
};


template<typename F1>
class Derivative<operations::Add<F1, functions::Const>> : public functions::Abstract
{
public:
	Derivative(const operations::Add<F1, functions::Const>& f)
		: m_df1(f.m_f1)
	{
	}

	double operator()(double x) override
	{
		double df1 = 0;

		if constexpr (std::is_pointer<Derivative<F1>>::value)
			df1 += (*m_df1)(x);
		else
			df1 += m_df1(x);

		return df1;
	}

	Derivative<F1> m_df1;

	typedef typename Derivative<F1>::Type Type;  // здесь мы оставляем только первый член выражения, потому что
	Type expression()                     // второй является константой
	{
		return m_df1.expression();
	}

};

template <typename F2>
class Derivative<operations::Add<functions::Const, F2>> : public functions::Abstract
{
public:
	Derivative(const operations::Add<functions::Const, F2>& f)
		: m_df2(f.m_f2)
	{
	}
	double operator()(double x) override
	{
		double df2 = 0;

		if constexpr (std::is_pointer<Derivative<F2>>::value)
			df2 += (*m_df2)(x);
		else
			df2 += m_df2(x);

		return df2;
	}

	Derivative<F2> m_df2;

	typedef typename Derivative<F2>::Type Type;   // Здесь второй
	Type expression()
	{
		return m_df2.expression();
	}
};

template<>
class Derivative<operations::Add<functions::Const, functions::Const>> : public functions::Abstract
{
public:
	Derivative(const operations::Add<functions::Const, functions::Const>&)
	{
	}

	double operator()(double) override
	{
		return 0;
	}

	typedef typename functions::Const Type;   // здесь две константы, ответ очевиден
	Type expression()
	{
		return functions::Const(0);
	}

};

template <>
class Derivative<operations::Add<functions::Simple, functions::Simple>> : public functions::Abstract
{
public:
	Derivative(const operations::Add<functions::Simple, functions::Simple>&)
	{
	}

	double operator()(double) override
	{
		return 2;
	}

	typedef functions::Const Type;  // Две функции "Simple" Ответ конечно же два!
	Type expression()
	{
		return functions::Const(2);
	}
};
