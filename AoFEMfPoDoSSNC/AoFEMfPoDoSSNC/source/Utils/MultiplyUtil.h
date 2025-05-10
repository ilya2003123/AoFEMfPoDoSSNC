#pragma once

#include "../Functions/functions.h"
#include "../Operations/operations.h"
#include <type_traits>

// Всё рассказал уже вроде

namespace utils
{
	template<typename F1, typename F2>
	operations::Multiply<F1, F2> operator*(const F1& f1, const F2& f2)
	{
		return operations::Multiply<F1, F2>(f1, f2);
	}

	template<typename F>
	typename std::enable_if<!std::is_arithmetic<F>::value, operations::Multiply<F, functions::Const>>::type
		operator*(const F& f, double value)
	{
		return operations::Multiply<F, functions::Const>(f, functions::Const(value));
	}

	template<typename F>
	typename std::enable_if<!std::is_arithmetic<F>::value, operations::Multiply<F, functions::Const>>::type
		operator*(const F& f, int value)
	{
		return operations::Multiply<F, functions::Const>(f, functions::Const(value));
	}

	template<typename F>
	typename std::enable_if<!std::is_arithmetic<F>::value, operations::Multiply<F, functions::Const>>::type
		operator*(double value, const F& f)
	{
		return operations::Multiply<F, functions::Const>(f, functions::Const(value));
	}

	template<typename F>
	typename std::enable_if<!std::is_arithmetic<F>::value, operations::Multiply<F, functions::Const>>::type
		operator*(int value, const F& f)
	{
		return operations::Multiply<F, functions::Const>(f, functions::Const(value));
	}

}