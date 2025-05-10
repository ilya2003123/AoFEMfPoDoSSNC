#pragma once

#include "../Functions/Abstract.h"

namespace operations  // Ещё немножко таких namespace-ов
{
	template<typename F1, typename F2>  // так же передаю в него 2 типа, вычитание же :D
	class Subtract : public functions::Abstract
	{
	public:
		typedef Subtract<F1, F2> Type;  // Обзываюсь
		Subtract(const F1& f1, const F2& f2)  // Конструктор
			:m_f1(f1), m_f2(f2)
		{
		}

		double operator()(double x) override  //  перегрузка
		{
			double f1 = 0;
			double f2 = 0;

			if constexpr (std::is_pointer<F1>::value)
				f1 += (*m_f1)(x);
			else
				f1 += m_f1(x);

			if constexpr (std::is_pointer_v<F2>)
				f2 += (*m_f2)(x);
			else
				f2 += m_f2(x);

			return f1 - f2;
		}

		F1 m_f1;  // Переменные
		F2 m_f2;

	};
}