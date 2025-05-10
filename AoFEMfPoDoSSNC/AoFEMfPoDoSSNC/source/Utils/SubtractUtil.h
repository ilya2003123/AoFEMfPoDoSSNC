#pragma once

#include "../Operations/operations.h"
#include "../Functions/functions.h"
#include <type_traits>

// Пройденный этап, глядите AddUtil.h

namespace utils
{
	template<typename F1, typename F2>
	operations::Subtract<F1, F2> operator-(const F1& f1, const F2& f2)
	{
		return operations::Subtract<F1, F2>(f1, f2);
	}

	template<typename F>
	typename std::enable_if<!std::is_arithmetic<F>::value, operations::Subtract<F, functions::Const>>::type
		operator-(const F& f, double value)
	{
		return operations::Subtract<F, functions::Const>(f, functions::Const(value));
	}

	template<typename F>
	typename std::enable_if<!std::is_arithmetic<F>::value, operations::Subtract<F, functions::Const>>::type
		operator-(const F& f, int value)
	{
		return operations::Subtract<F, functions::Const>(f, functions::Const(value));
	}

	template<typename F>
	typename std::enable_if<!std::is_arithmetic<F>::value, operations::Subtract<functions::Const, F>>::type
		operator-(double value, const F& f)
	{
		return operations::Subtract<functions::Const, F>(functions::Const(value), f);
	}

	template<typename F>
	typename std::enable_if<!std::is_arithmetic<F>::value, operations::Subtract<functions::Const, F>>::type
		operator-(int value, const F& f)
	{
		return operations::Subtract<functions::Const, F>(functions::Const(value), f);
	}

}