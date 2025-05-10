#pragma once

#include "../Functions/Abstract.h"

namespace operations  // Ну тут уже немножко другой namespace, нововведения произошли очень быстро
{
	template<typename F1, typename F2>  // тут передаём уже 2 класса! ну или же, по правильному сказать 2 типа
	class Add : public functions::Abstract // Класс суммы
	{
	public:  // обзываю тот тип, который получится из сложения двух типов, которые я передам в этот класс
		typedef Add<F1, F2> Type;
		Add(const F1& f1, const F2& f2) 
			:m_f1(f1), m_f2(f2)
		{
		}

		double operator()(double x) override  // Ну, с такими выражениями мы уже знакомы
		{
			double f1 = 0;
			double f2 = 0;
			if constexpr (std::is_pointer_v<F1>)
				f1 += (*m_f1)(x);
			else
				f1 += m_f1(x);

			if constexpr (std::is_pointer_v<F2>)
				f2 += (*m_f2)(x);
			else
				f2 += m_f2(x);

			return f1 + f2;
		}

		F1 m_f1; // Одна переменная одного типа, вторая другого
		F2 m_f2;

	};
}