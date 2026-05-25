#pragma once

#include <memory>
#include <type_traits>

namespace functions
{
	class Abstract
	{
	public:
		virtual ~Abstract() = default;
		virtual double operator()(double x) = 0;
	};

	template<typename F>
	double evaluate(F& f, double x)
	{
		if constexpr (std::is_pointer_v<F>)
		{
			return (*f)(x);
		}
		else
		{
			return f(x);
		}
	}

	template<typename T>
	double evaluate(std::shared_ptr<T>& f, double x)
	{
		return (*f)(x);
	}

	template<typename T>
	double evaluate(const std::shared_ptr<T>& f, double x)
	{
		return (*f)(x);
	}
}
