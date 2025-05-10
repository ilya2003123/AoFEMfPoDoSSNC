#pragma once

#include "../Functions/functions.h"
#include "../Operations/operations.h"

namespace utils
{
	template<typename F>
	operations::UnarMinus<F> operator-(const F& f)
	{
		return operations::UnarMinus<F>(f);
	}
}