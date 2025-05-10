#pragma once

#include <cmath>
#include "Abstract.h"

namespace functions  // тот же самый namespace
{
	template<typename F> // шаблончик, вот тут то и начинается веселье, передаю в него тип F (может быть любой)
	class Power : public functions::Abstract
	{
	public:   // обязательно всё public. Ошибки сыпятся!
		typedef Power<F> Type;  // Обзываю свою степень в которую передал какой-то тип Type-ом, тоже тёска
		Power(const F& f, const double n)   // Продолжаю эпоху минимализма
			: m_f(f), m_n(n)
		{
		}

		double operator()(double x) override  // С этим уже вроде знакомы
		{
			if constexpr (std::is_pointer<F>::value)
				return pow((*m_f)(x), m_n);
			else
				return pow(m_f(x), m_n);
		}

		F m_f;        // Одна переменная типа F
		double m_n;   // Вторая double

	};
}