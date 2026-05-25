#pragma once

#include <functional>

// Numerical derivative for arbitrary callable expression nodes.
template<typename F>
class Derivative : public functions::Abstract
{
public:
	Derivative(const F& f, double dx = 5e-6)
		: m_f(f), m_dx(dx)
	{
	}

	double operator()(double x) override
	{
		const double fx1 = functions::evaluate(m_f, x + m_dx);
		const double fx2 = functions::evaluate(m_f, x);
		return (fx1 - fx2) / m_dx;
	}

	using Type = std::function<double(double)>;

	Type expression()
	{
		return [this](double x)
		{
			return (functions::evaluate(m_f, x + m_dx) - functions::evaluate(m_f, x)) / m_dx;
		};
	}

private:
	F m_f;
	double m_dx;
};
