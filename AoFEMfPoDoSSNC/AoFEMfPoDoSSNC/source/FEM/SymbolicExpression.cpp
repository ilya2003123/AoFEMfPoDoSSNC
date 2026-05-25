#include "SymbolicExpression.h"

#include <algorithm>
#include <cmath>
#include <cctype>
#include <cstdlib>
#include <iomanip>
#include <map>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace symbolic
{
	std::string formatNumber(double value)
	{
		if (std::abs(value) < 1e-12)
		{
			value = 0.0;
		}

		std::ostringstream out;
		out << std::setprecision(12) << value;
		return out.str();
	}

	SymbolicExpression number(double value)
	{
		SymbolicExpression expression;
		expression.kind = SymbolicExpression::Kind::Number;
		expression.number = value;
		expression.token = formatNumber(value);
		return expression;
	}

	SymbolicExpression variable()
	{
		SymbolicExpression expression;
		expression.kind = SymbolicExpression::Kind::Variable;
		expression.token = "x";
		return expression;
	}

	bool isNumber(const SymbolicExpression& expression, double value)
	{
		return expression.kind == SymbolicExpression::Kind::Number
			&& std::abs(expression.number - value) < 1e-12;
	}

	SymbolicExpression unary(std::string token, SymbolicExpression arg)
	{
		if (token == "-")
		{
			if (arg.kind == SymbolicExpression::Kind::Number)
			{
				return number(-arg.number);
			}
			if (arg.kind == SymbolicExpression::Kind::Unary && arg.token == "-")
			{
				return arg.args[0];
			}
		}

		SymbolicExpression expression;
		expression.kind = SymbolicExpression::Kind::Unary;
		expression.token = std::move(token);
		expression.args.push_back(std::move(arg));
		return expression;
	}

	SymbolicExpression binary(std::string token, SymbolicExpression left, SymbolicExpression right)
	{
		if (left.kind == SymbolicExpression::Kind::Number && right.kind == SymbolicExpression::Kind::Number)
		{
			if (token == "+") return number(left.number + right.number);
			if (token == "-") return number(left.number - right.number);
			if (token == "*") return number(left.number * right.number);
			if (token == "/" && std::abs(right.number) > 1e-12) return number(left.number / right.number);
			if (token == "^") return number(std::pow(left.number, right.number));
		}

		if (token == "+")
		{
			if (isNumber(left, 0.0)) return right;
			if (isNumber(right, 0.0)) return left;
		}
		else if (token == "-")
		{
			if (isNumber(right, 0.0)) return left;
			if (isNumber(left, 0.0)) return unary("-", right);
		}
		else if (token == "*")
		{
			if (isNumber(left, 0.0) || isNumber(right, 0.0)) return number(0.0);
			if (isNumber(left, 1.0)) return right;
			if (isNumber(right, 1.0)) return left;
			if (isNumber(left, -1.0)) return unary("-", right);
			if (isNumber(right, -1.0)) return unary("-", left);

			if (left.kind == SymbolicExpression::Kind::Number
				&& right.kind == SymbolicExpression::Kind::Binary
				&& right.token == "*"
				&& right.args[0].kind == SymbolicExpression::Kind::Number)
			{
				return binary("*", number(left.number * right.args[0].number), right.args[1]);
			}
			if (right.kind == SymbolicExpression::Kind::Number
				&& left.kind == SymbolicExpression::Kind::Binary
				&& left.token == "*"
				&& left.args[0].kind == SymbolicExpression::Kind::Number)
			{
				return binary("*", number(right.number * left.args[0].number), left.args[1]);
			}
		}
		else if (token == "/")
		{
			if (isNumber(left, 0.0)) return number(0.0);
			if (isNumber(right, 1.0)) return left;
		}
		else if (token == "^")
		{
			if (isNumber(right, 0.0)) return number(1.0);
			if (isNumber(right, 1.0)) return left;
		}

		SymbolicExpression expression;
		expression.kind = SymbolicExpression::Kind::Binary;
		expression.token = std::move(token);
		expression.args.push_back(std::move(left));
		expression.args.push_back(std::move(right));
		return expression;
	}

	int precedence(const SymbolicExpression& expression)
	{
		if (expression.kind == SymbolicExpression::Kind::Binary)
		{
			if (expression.token == "+" || expression.token == "-") return 1;
			if (expression.token == "*" || expression.token == "/") return 2;
			if (expression.token == "^") return 4;
		}
		if (expression.kind == SymbolicExpression::Kind::Unary) return 3;
		return 5;
	}

	std::string toString(const SymbolicExpression& expression, int parentPrecedence)
	{
		std::string result;
		const int ownPrecedence = precedence(expression);

		switch (expression.kind)
		{
		case SymbolicExpression::Kind::Number:
		case SymbolicExpression::Kind::Variable:
			result = expression.token;
			break;
		case SymbolicExpression::Kind::Unary:
			if (expression.token == "-")
			{
				result = "-" + toString(expression.args[0], ownPrecedence);
			}
			else
			{
				result = expression.token + "(" + toString(expression.args[0]) + ")";
			}
			break;
		case SymbolicExpression::Kind::Binary:
		{
			const bool parenthesizeRight =
				(expression.token == "-" || expression.token == "/")
				&& expression.args[1].kind == SymbolicExpression::Kind::Binary
				&& precedence(expression.args[1]) <= ownPrecedence;
			const bool parenthesizeLeft =
				expression.token == "^"
				&& expression.args[0].kind == SymbolicExpression::Kind::Binary
				&& expression.args[0].token == "^";

			const std::string left = parenthesizeLeft
				? "(" + toString(expression.args[0]) + ")"
				: toString(expression.args[0], ownPrecedence);
			const std::string right = parenthesizeRight
				? "(" + toString(expression.args[1]) + ")"
				: toString(expression.args[1], ownPrecedence + (expression.token == "^" ? -1 : 0));

			result = left + expression.token + right;
			break;
		}
		case SymbolicExpression::Kind::Piecewise:
			result = "pw(";
			for (int i = 0; i < static_cast<int>(expression.args.size()); ++i)
			{
				if (i > 0)
				{
					result += ", ";
				}
				result += toString(expression.args[i]);
			}
			result += ")";
			break;
		}

		if (ownPrecedence < parentPrecedence)
		{
			return "(" + result + ")";
		}
		return result;
	}

	SymbolicExpression derivative(const SymbolicExpression& expression)
	{
		switch (expression.kind)
		{
		case SymbolicExpression::Kind::Number:
			return number(0.0);
		case SymbolicExpression::Kind::Variable:
			return number(1.0);
		case SymbolicExpression::Kind::Piecewise:
		{
			SymbolicExpression result;
			result.kind = SymbolicExpression::Kind::Piecewise;
			result.token = "pw";
			for (int i = 0; i < static_cast<int>(expression.args.size()); ++i)
			{
				const bool isBorder = (i % 2 == 0) && (i + 1 < static_cast<int>(expression.args.size()));
				result.args.push_back(isBorder ? expression.args[i] : derivative(expression.args[i]));
			}
			return result;
		}
		case SymbolicExpression::Kind::Unary:
		{
			const SymbolicExpression& arg = expression.args[0];
			SymbolicExpression dArg = derivative(arg);

			if (expression.token == "-") return unary("-", dArg);
			if (expression.token == "sin") return binary("*", unary("cos", arg), dArg);
			if (expression.token == "cos") return unary("-", binary("*", unary("sin", arg), dArg));
			if (expression.token == "tg" || expression.token == "tan")
			{
				return binary("/", dArg, binary("^", unary("cos", arg), number(2.0)));
			}
			if (expression.token == "ctg" || expression.token == "cot")
			{
				return unary("-", binary("/", dArg, binary("^", unary("sin", arg), number(2.0))));
			}
			if (expression.token == "sqrt")
			{
				return binary("/", dArg, binary("*", number(2.0), unary("sqrt", arg)));
			}
			if (expression.token == "sqr")
			{
				return binary("*", binary("*", number(2.0), arg), dArg);
			}
			if (expression.token == "exp")
			{
				return binary("*", unary("exp", arg), dArg);
			}
			if (expression.token == "log")
			{
				return binary("/", dArg, arg);
			}
			throw std::runtime_error("Unsupported function for generated f display: " + expression.token);
		}
		case SymbolicExpression::Kind::Binary:
		{
			const SymbolicExpression& left = expression.args[0];
			const SymbolicExpression& right = expression.args[1];
			SymbolicExpression dLeft = derivative(left);
			SymbolicExpression dRight = derivative(right);

			if (expression.token == "+") return binary("+", dLeft, dRight);
			if (expression.token == "-") return binary("-", dLeft, dRight);
			if (expression.token == "*")
			{
				return binary("+", binary("*", dLeft, right), binary("*", left, dRight));
			}
			if (expression.token == "/")
			{
				return binary(
					"/",
					binary("-", binary("*", dLeft, right), binary("*", left, dRight)),
					binary("^", right, number(2.0)));
			}
			if (expression.token == "^")
			{
				if (right.kind == SymbolicExpression::Kind::Number)
				{
					return binary(
						"*",
						binary("*", number(right.number), binary("^", left, number(right.number - 1.0))),
						dLeft);
				}
				return binary(
					"*",
					expression,
					binary(
						"+",
						binary("*", dRight, unary("log", left)),
						binary("/", binary("*", right, dLeft), left)));
			}
			throw std::runtime_error("Unsupported operator for generated f display: " + expression.token);
		}
		}

		return number(0.0);
	}
}

