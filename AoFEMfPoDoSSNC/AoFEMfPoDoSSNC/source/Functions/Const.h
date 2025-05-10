#pragma once

#include "Abstract.h"

namespace functions   // Свои функции желательно облеплять в namespace, вдруг у какой-то функции
{					  // будет такое же название O_o
	class Const : public functions::Abstract
	{
	public:
		typedef Const Type;  // Обзываю Const Type-ом, понадобится дальше
		Const(const double x)  // Обычный минималистичный конструктор, для типа double
			: m_const(x)
		{
		}

		Const(const int x)   // Тоже конструктор, только целочисленный
			: m_const(x)
		{
		}

		double operator()(double) override  // перегружаю оператор вызова функции
		{                                // возвращать будет double и передаваться в неё будет тоже double
			return m_const;
		}                                // Можно написать double operator()(double x) но зачем ? :D
		                                 // x я там не использую

		operator double()  // Очень хитрая штука, говорю что мой класс может быть double (оператор преобразования)
		{                  // Преобразую свой класс в double
			return m_const;
		}

		double m_const;   // Минималистичная переменная, отвечающая за константу
	};

}