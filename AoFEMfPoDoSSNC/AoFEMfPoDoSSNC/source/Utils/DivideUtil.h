#pragma once

#include "../Functions/functions.h"
#include "../Operations/operations.h"
#include <type_traits>

namespace utils
{
	template<typename F1, typename F2>
	operations::Divide<F1, F2> operator/(const F1& f1, const F2& f2)
	{
		return operations::Divide<F1, F2>(f1, f2);
	}

	template<typename F>
	typename std::enable_if<!std::is_arithmetic<F>::value, operations::Divide<F, functions::Const>>::type
		operator/(const F& f, double value)
	{
		return operations::Divide<F, functions::Const>(f, functions::Const(value));
	}

	template<typename F>
	typename std::enable_if<!std::is_arithmetic<F>::value, operations::Divide<F, functions::Const>>::type
		operator/(const F& f, int value)
	{
		return operations::Divide<F, functions::Const>(f, functions::Const(value));
	}

	template<typename F>
	typename std::enable_if<!std::is_arithmetic<F>::value, operations::Divide<functions::Const, F>>::type
		operator/(double value, const F& f)
	{
		return operations::Divide<functions::Const, F>(functions::Const(value), f);
	}

	template<typename F>
	typename std::enable_if<!std::is_arithmetic<F>::value, operations::Divide<functions::Const, F>>::type
		operator/(int value, const F& f)
	{
		return operations::Divide<functions::Const, F>(functions::Const(value), f);
	}

}