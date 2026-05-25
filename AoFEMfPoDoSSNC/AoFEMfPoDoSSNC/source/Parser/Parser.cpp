#include "Parser.h"

#include "../Derivative/derivative.h"
#include "../Operations/operations.h"
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <stdexcept>
#include <utility>

namespace
{
	using FunctionPtr = std::shared_ptr<functions::Abstract>;
}

Expression::Expression(std::string token, std::shared_ptr<functions::Abstract> f)
	: token(std::move(token)), func(std::move(f))
{
}

Parser::Parser(const char* input)
	: input(input), begin(input)
{
}

void Parser::skipSpaces()
{
	while (std::isspace(static_cast<unsigned char>(*input)))
	{
		++input;
	}
}

bool Parser::match(char expected)
{
	skipSpaces();
	return *input == expected;
}

bool Parser::consume(char expected)
{
	if (!match(expected))
	{
		return false;
	}
	++input;
	return true;
}

std::string Parser::parseIdentifier()
{
	skipSpaces();
	std::string result;
	while (std::isalpha(static_cast<unsigned char>(*input)) || *input == '_')
	{
		result.push_back(*input++);
	}
	return result;
}

Expression Parser::makeConstant(double value)
{
	auto node = std::make_shared<functions::Const>(value);
	Expression result("number", node);
	result.number = value;
	result.isNumber = true;
	return result;
}

Expression Parser::parseNumber()
{
	skipSpaces();
	char* end = nullptr;
	double value = std::strtod(input, &end);
	if (end == input)
	{
		throw std::runtime_error("Expected number");
	}
	input = end;
	return makeConstant(value);
}

Expression Parser::makeUnary(const std::string& token, const Expression& arg)
{
	Expression result;
	result.token = token;

	if (token == "+")
	{
		return arg;
	}
	if (token == "-")
	{
		result.func = std::make_shared<operations::UnarMinus<FunctionPtr>>(arg.func);
	}
	else if (token == "sin")
	{
		result.func = std::make_shared<functions::Sin<FunctionPtr>>(arg.func);
	}
	else if (token == "cos")
	{
		result.func = std::make_shared<functions::Cos<FunctionPtr>>(arg.func);
	}
	else if (token == "tg" || token == "tan")
	{
		result.func = std::make_shared<functions::Tan<FunctionPtr>>(arg.func);
	}
	else if (token == "ctg" || token == "cot")
	{
		result.func = std::make_shared<functions::Cot<FunctionPtr>>(arg.func);
	}
	else if (token == "asin" || token == "arcsin")
	{
		result.func = std::make_shared<functions::ArcSin<FunctionPtr>>(arg.func);
	}
	else if (token == "acos" || token == "arccos")
	{
		result.func = std::make_shared<functions::ArcCos<FunctionPtr>>(arg.func);
	}
	else if (token == "atg" || token == "atan" || token == "arctan")
	{
		result.func = std::make_shared<functions::ArcTan<FunctionPtr>>(arg.func);
	}
	else if (token == "actg" || token == "acot" || token == "arccot")
	{
		result.func = std::make_shared<functions::ArcCot<FunctionPtr>>(arg.func);
	}
	else if (token == "sqrt")
	{
		result.func = std::make_shared<functions::Power<FunctionPtr>>(arg.func, 0.5);
	}
	else if (token == "sqr")
	{
		result.func = std::make_shared<functions::Power<FunctionPtr>>(arg.func, 2.0);
	}
	else if (token == "dx")
	{
		result.func = std::make_shared<Derivative<FunctionPtr>>(arg.func);
	}
	else
	{
		throw std::runtime_error("Unknown unary function: " + token);
	}

	return result;
}

Expression Parser::makeBinary(const std::string& token, const Expression& left, const Expression& right)
{
	Expression result;
	result.token = token;

	if (token == "+")
	{
		result.func = std::make_shared<operations::Add<FunctionPtr, FunctionPtr>>(left.func, right.func);
	}
	else if (token == "-")
	{
		result.func = std::make_shared<operations::Subtract<FunctionPtr, FunctionPtr>>(left.func, right.func);
	}
	else if (token == "*")
	{
		result.func = std::make_shared<operations::Multiply<FunctionPtr, FunctionPtr>>(left.func, right.func);
	}
	else if (token == "/")
	{
		result.func = std::make_shared<operations::Divide<FunctionPtr, FunctionPtr>>(left.func, right.func);
	}
	else if (token == "^")
	{
		result.func = std::make_shared<functions::Exponent_Power<FunctionPtr, FunctionPtr>>(left.func, right.func);
	}
	else
	{
		throw std::runtime_error("Unknown binary operator: " + token);
	}

	return result;
}

std::vector<Expression> Parser::parseArguments()
{
	std::vector<Expression> args;
	if (consume(')'))
	{
		return args;
	}

	for (;;)
	{
		args.push_back(parseAddSub());
		if (consume(')'))
		{
			return args;
		}
		if (!consume(','))
		{
			throw std::runtime_error("Expected ',' or ')'");
		}
	}
}

Expression Parser::parsePrimary()
{
	skipSpaces();
	char current = *input;

	if (current == '\0')
	{
		throw std::runtime_error("Unexpected end of expression");
	}

	if (std::isdigit(static_cast<unsigned char>(current)) || current == '.')
	{
		return parseNumber();
	}

	if (consume('('))
	{
		Expression inside = parseAddSub();
		if (!consume(')'))
		{
			throw std::runtime_error("Expected ')'");
		}
		return inside;
	}

	std::string id = parseIdentifier();
	if (id.empty())
	{
		throw std::runtime_error("Unexpected character in expression");
	}

	if (id == "x")
	{
		return Expression("x", std::make_shared<functions::Simple>());
	}
	if (id == "pi")
	{
		return makeConstant(std::acos(-1.0));
	}
	if (id == "e")
	{
		return makeConstant(std::exp(1.0));
	}

	std::vector<Expression> args;
	if (consume('('))
	{
		args = parseArguments();
	}
	else
	{
		args.push_back(parseUnary());
	}

	if (id == "pow")
	{
		if (args.size() != 2)
		{
			throw std::runtime_error("pow expects 2 arguments");
		}
		if (args[1].isNumber)
		{
			Expression result("pow", std::make_shared<functions::Power<FunctionPtr>>(args[0].func, args[1].number));
			return result;
		}
		return makeBinary("^", args[0], args[1]);
	}

	if (id == "pw" || id == "piecewise")
	{
		if (args.size() < 3 || args.size() % 2 == 0)
		{
			throw std::runtime_error("pw expects: pw(border1, expr1, ..., lastExpr)");
		}

		std::vector<double> borders;
		std::vector<FunctionPtr> pieces;
		for (int i = 0; i < static_cast<int>(args.size()) - 1; i += 2)
		{
			if (!args[i].isNumber)
			{
				throw std::runtime_error("piecewise border must be a number");
			}
			borders.push_back(args[i].number);
			pieces.push_back(args[i + 1].func);
		}
		pieces.push_back(args.back().func);

		return Expression(id, std::make_shared<functions::Piecewise>(borders, pieces));
	}

	if (id == "exp")
	{
		if (args.size() == 1)
		{
			return Expression("exp", std::make_shared<functions::Exponent<FunctionPtr>>(std::exp(1.0), args[0].func));
		}
		if (args.size() == 2 && args[0].isNumber)
		{
			return Expression("exp", std::make_shared<functions::Exponent<FunctionPtr>>(args[0].number, args[1].func));
		}
		throw std::runtime_error("exp expects exp(expr) or exp(number, expr)");
	}

	if (id == "log")
	{
		if (args.size() == 1)
		{
			auto base = std::make_shared<functions::Const>(std::exp(1.0));
			return Expression("log", std::make_shared<functions::Logarithm<FunctionPtr, FunctionPtr>>(base, args[0].func));
		}
		if (args.size() == 2)
		{
			return Expression("log", std::make_shared<functions::Logarithm<FunctionPtr, FunctionPtr>>(args[0].func, args[1].func));
		}
		throw std::runtime_error("log expects 1 or 2 arguments");
	}

	if (args.size() != 1)
	{
		throw std::runtime_error(id + " expects 1 argument");
	}
	return makeUnary(id, args[0]);
}

Expression Parser::parseUnary()
{
	if (consume('+'))
	{
		return parseUnary();
	}
	if (consume('-'))
	{
		return makeUnary("-", parseUnary());
	}
	return parsePower();
}

Expression Parser::parsePower()
{
	Expression left = parsePrimary();
	if (consume('^'))
	{
		Expression right = parseUnary();
		return makeBinary("^", left, right);
	}
	return left;
}

Expression Parser::parseMulDiv()
{
	Expression left = parseUnary();
	for (;;)
	{
		if (consume('*'))
		{
			left = makeBinary("*", left, parseUnary());
		}
		else if (consume('/'))
		{
			left = makeBinary("/", left, parseUnary());
		}
		else
		{
			return left;
		}
	}
}

Expression Parser::parseAddSub()
{
	Expression left = parseMulDiv();
	for (;;)
	{
		if (consume('+'))
		{
			left = makeBinary("+", left, parseMulDiv());
		}
		else if (consume('-'))
		{
			left = makeBinary("-", left, parseMulDiv());
		}
		else
		{
			return left;
		}
	}
}

Expression Parser::parse()
{
	begin = input;
	Expression result = parseAddSub();
	skipSpaces();
	if (*input != '\0')
	{
		throw std::runtime_error("Unexpected tail near position " + std::to_string(input - begin));
	}
	return result;
}

double eval(const Expression& e, double x)
{
	return (*e.func)(x);
}
