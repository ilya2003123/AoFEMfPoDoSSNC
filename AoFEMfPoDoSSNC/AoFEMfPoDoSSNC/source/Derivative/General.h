#pragma once

#include <functional>  // functional для std::function<>

// тут разворачивается неописуемая магия, нет include "function.h", но он всё равно видит Abstract

// Здесь решил почему-то обойтись без namespace, Derivative я использовать не буду чтоб прям часто, поэтому
// был сделан такой рискованный шаг, если захотите, могу обернуть в namespace

template<typename F> // Это так сказать, БАЗОВЫЙ класс производных
class Derivative : public functions::Abstract
{
public:
	Derivative(const F& f, double dx = 5e-6) // конструктор, дельту беру за (1 / 5000000)
		: m_f(f), m_dx(dx)
	{
	}

	double operator()(double x) override
	{
		double fx1 = 0;
		double fx2 = 0;

		if constexpr (std::is_pointer<F>::value)
			fx1 += (*m_f)(x + m_dx);
		else
			fx1 += m_f(x + m_dx);

		if constexpr (std::is_pointer<F>::value)
			fx2 += (*m_f)(x);
		else
			fx2 += m_f(x);

		return ((fx1 - fx2) / m_dx);  // численное вычисление производных
	}

	F m_f;
	double m_dx;

	// std::function это шаблонный класс, который является "оболочкой" (контейнером) для каких-то вещей 
	typedef std::function<double(double)> Type; // (функции, лямба-выражения и т.д.) А мы же все функции перегружали
	Type expression()  // в double (double) поэтому, грубо говоря, она хватает сюда функцию, любую и обзывает её
	{
		return [&](double x)  // лямбда-выражение [](){}, [&] - может хватать всё что выше неё, по ссылке, 
		{          // чтобы потом использовать, здесь оно надо чтобы заработала передача нашего икса
			return ((m_f(x + m_dx) - m_f(x)) / m_dx);  // также тут return функции, и сразу её реализация
		};							// Опережаю вопрос: "Почему реализацию прям тут?" - Таков синтаксис :D
	}								// Красота, не правда ли?
};