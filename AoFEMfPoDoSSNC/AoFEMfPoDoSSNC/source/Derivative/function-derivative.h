#pragma once

#include "General.h"

// Эта шаблонная функция создана только ради красоты, чтобы можно было вызывать её, а не мучаться с expression()
// Поверьте мне, лучше так, чем сидеть и думать, почему проект собирается и даже компилируется, а потом 
// вылетает ошибка от того, что шаг дельта(m_dx) равен 6.7923578274e+61

template <typename F>
typename /*Derivative<F>::Type*/ auto derivative(F f)
{
	return Derivative<F>(f).expression();
}