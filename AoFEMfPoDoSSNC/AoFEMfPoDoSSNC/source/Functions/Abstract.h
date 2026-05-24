#pragma once

#include <type_traits>

namespace functions
{
	class Abstract
	{
	public:
		virtual ~Abstract() = default;
		virtual double operator()(double x) = 0;
	};
}
