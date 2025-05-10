#pragma once

#include"../Functions/functions.h"
#include "General.h"

template<>  // Пустые <> "скобочки" говорят нам о частичной специализации этого шаблона
class Derivative<functions::Simple> : public functions::Abstract  
{   // Грубо говоря, мы перегружаем класс наших производных, для Simple
public:
	typedef functions::Const Type;  // Производная понятное дело будет константа же o_O
	Derivative(functions::Simple) {}  // Вот такой минималистичный конструктор

	double operator()(double) override // знакомые вещи
	{
		return 1;
	}

	Type expression()
	{
		return functions::Const(1);   // тип возвращаемый будет Const и возвращаем мы Const, всё честно
	}

};