#pragma once

#include "../Operations/operations.h"
#include "../Functions/functions.h"
#include <type_traits>  // для enable_if, а также is_arithmetic

// тут всё очень страшно написано, но всё очень просто

namespace utils  // ещё один блок, где будут свои функции
{	
	template<typename F1, typename F2>  // Передаю 2 типа
	operations::Add<F1, F2> operator+(const F1& f1, const F2& f2) // перегружаю оператор + для 2 каких-то типов
	{
		return operations::Add<F1, F2>(f1, f2);  // Возвращаю просто класс
	}

// здесь предстоит сложить константу и какой-то тип, сложение может быть 4 типов (double + F, F + double)
// (int + F, F + int) перегружаем, НО перед этим мы спросим (enable_if) является ли наш класс F
// НЕ численным типом, F в лёгкую может быть константой, "::value" -> возвращает true или false,
// После чего мы закидываем туда операцию сложения типа F и Const, ну и потом спрашиваем тип этого всего 
// дела "::type", чтобы было понятно какой тип возвращать. String сюда передать нельзя, проверенно.

	template<typename F> // передаю 1 тип
	typename std::enable_if<!std::is_arithmetic<F>::value, operations::Add<F, functions::Const>>::type
		operator+(double value, const F& f)
	{
		return operations::Add<F, functions::Const>(f, functions::Const(value));
	}

	template<typename F>
	typename std::enable_if<!std::is_arithmetic<F>::value, operations::Add<F, functions::Const>>::type
		operator+(const F& f, double value)
	{
		return operations::Add<F, functions::Const>(f, functions::Const(value));
	}

	template<typename F>
	typename std::enable_if<!std::is_arithmetic<F>::value, operations::Add<F, functions::Const>>::type
		operator+(int value, const F& f)
	{
		return operations::Add<F, functions::Const>(f, functions::Const(value));
	}

	template<typename F>
	typename std::enable_if<!std::is_arithmetic<F>::value, operations::Add<F, functions::Const>>::type
		operator+(const F& f, int value)
	{
		return operations::Add<F, functions::Const>(f, functions::Const(value));
	}

}