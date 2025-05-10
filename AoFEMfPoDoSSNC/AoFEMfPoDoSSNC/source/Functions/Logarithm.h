#pragma once

#include "Abstract.h"
#include <cmath>	// Думаю, пора рассказать про прекрасную стандартную библиотеку языка С++
					// Там есть нужные нам функции, всё (-_-)
	

namespace functions 
{
	template<typename F1, typename F2>   // Тут всё ясненько вроде
	class Logarithm : public functions::Abstract
	{
	public:
		typedef Logarithm<F1, F2> Type;  // Тут вроде тоже
		Logarithm(const F1& f1, const F2& f2) // Вот тут могут возникнуть вопросики, я предвидел это, поэтому
			:m_f1(f1), m_f2(f2) // сразу говорю. у нас в "cmath" есть только
		{    // натуральный логарифм, нам такое не подходит, мы же хотим ещё какой-то, для нас существует много
		}	 // логарифмов, поэтому мы по одному из свойств логарифма можем получить любой, нам надо разделить
	    // два логарифма с одинаковыми основаниями, короче вы и так это знаете лучше меня. Точка после 1 нужна!	

		double operator()(double x) override
		{
			double factor = 0;
			double fx = 0;
			if constexpr (std::is_pointer<F2>::value)
				fx = log((*m_f2)(x));
			else
				fx = log(m_f2(x)); 

			if constexpr (std::is_pointer<F1>::value)
				factor = (1. / log((*m_f1)(x)));
			else
				factor = (1. / log(m_f1(x)));

			return fx * factor;
			// Вот тут мы это и делаем, внимательно посмотрите чему равны переменные
		}

		F1 m_f1;
		F2 m_f2;
	};
}