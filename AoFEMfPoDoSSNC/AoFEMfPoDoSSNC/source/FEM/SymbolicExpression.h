#pragma once

#include <string>
#include <vector>

namespace symbolic
{
	struct SymbolicExpression
	{
		enum class Kind
		{
			Number,
			Variable,
			Unary,
			Binary,
			Piecewise
		};

		Kind kind = Kind::Number;
		std::string token;
		double number = 0.0;
		std::vector<SymbolicExpression> args;
	};

	SymbolicExpression parseExpression(const char* input);
	std::string formatNumber(double value);
	SymbolicExpression number(double value);
	SymbolicExpression variable();
	SymbolicExpression unary(std::string token, SymbolicExpression arg);
	SymbolicExpression binary(std::string token, SymbolicExpression left, SymbolicExpression right);
	SymbolicExpression derivative(const SymbolicExpression& expression);
	std::string toString(const SymbolicExpression& expression, int parentPrecedence = 0);
	std::string toSimplifiedString(const SymbolicExpression& expression);
}
