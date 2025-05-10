#pragma once

#include <cmath>
#include "Abstract.h"

namespace functions // Ну тут уж слов нет
{
	template<typename F>  // Вроде тоже всё ясно
	class Exponent : public functions::Abstract
	{
	public:  // Закрепим, обзываю свой класс показательной функции в которую передаю тип F, Type-ом
		typedef Exponent<F> Type;  
		Exponent(const double base, const F& f)
			:m_base(base), m_f(f)
		{
		}

		double operator()(double x) override  // Ну тут говорить нечего, я считаю
		{
			if constexpr (std::is_pointer<F>::value)
				return pow(m_base, (*m_f)(x));
			else
				return pow(m_base, m_f(x));   // используем стандартные формулы языка,
		}								  // мы што умные штоли с нуля ВСЁ писать :D

		double m_base;  // здесь основание double
		F m_f;          // А вот здесь будет наш тип F (предположительно какая-то функция)

	};
}