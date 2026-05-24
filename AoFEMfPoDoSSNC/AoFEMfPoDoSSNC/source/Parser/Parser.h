#pragma once

#pragma warning(disable: 26495)

#include "../Functions/functions.h"
#include "../Operations/operations.h"
#include "../Derivative/derivative.h"
#include <cmath>
#include <string>
#include <vector>
#include <cctype>
#include <cstring>
#include <cstdlib>
#include <memory>
#include <stdexcept>
#include <utility>

inline double inputx = 0;

namespace functions
{
	class PiecewiseRuntime : public functions::Abstract
	{
	public:
		PiecewiseRuntime(std::vector<double> borders, std::vector<std::shared_ptr<functions::Abstract>> pieces)
			: m_borders(std::move(borders)), m_pieces(std::move(pieces))
		{
		}

		double operator()(double x) override
		{
			for (int i = 0; i < static_cast<int>(m_borders.size()); ++i)
			{
				if (x <= m_borders[i] + 1e-12)
				{
					return (*m_pieces[i])(x);
				}
			}
			return (*m_pieces.back())(x);
		}

	private:
		std::vector<double> m_borders;
		std::vector<std::shared_ptr<functions::Abstract>> m_pieces;
	};
}

struct Expression 
{
	Expression() = default;
	Expression(std::string token, std::shared_ptr<functions::Abstract> f)
		: token(std::move(token)), func(std::move(f))
	{
	}

	std::string token;
	std::shared_ptr<functions::Abstract> func;
	std::vector<std::shared_ptr<functions::Abstract>> keepAlive;
	double number = 0;
	bool isNumber = false;
};

class Parser 
{
public:
	explicit Parser(const char* input)
		: input(input), begin(input)
	{
	}
	Expression parse();

private:
	void skipSpaces();
	bool match(char expected);
	bool consume(char expected);
	char peek();
	std::string parseIdentifier();
	Expression parseNumber();
	Expression parsePrimary();
	Expression parseUnary();
	Expression parsePower();
	Expression parseMulDiv();
	Expression parseAddSub();
	std::vector<Expression> parseArguments();
	Expression makeConstant(double value);
	Expression makeUnary(const std::string& token, const Expression& arg);
	Expression makeBinary(const std::string& token, const Expression& left, const Expression& right);
	void keep(Expression& result, const Expression& child);
	void keep(Expression& result, const std::shared_ptr<functions::Abstract>& ptr);

	const char* input;
	const char* begin;
};

void Parser::skipSpaces()
{
	while (std::isspace(static_cast<unsigned char>(*input))) ++input;
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

char Parser::peek()
{
	skipSpaces();
	return *input;
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

void Parser::keep(Expression& result, const std::shared_ptr<functions::Abstract>& ptr)
{
	if (ptr)
	{
		result.keepAlive.push_back(ptr);
	}
}

void Parser::keep(Expression& result, const Expression& child)
{
	keep(result, child.func);
	result.keepAlive.insert(result.keepAlive.end(), child.keepAlive.begin(), child.keepAlive.end());
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
		result.func = std::make_shared<operations::UnarMinus<functions::Abstract*>>(arg.func.get());
	}
	else if (token == "sin")
	{
		result.func = std::make_shared<functions::Sinus<functions::Abstract*>>(arg.func.get());
	}
	else if (token == "cos")
	{
		result.func = std::make_shared<functions::Cosinus<functions::Abstract*>>(arg.func.get());
	}
	else if (token == "tg")
	{
		result.func = std::make_shared<functions::Tangent<functions::Abstract*>>(arg.func.get());
	}
	else if (token == "ctg")
	{
		result.func = std::make_shared<functions::Cotangent<functions::Abstract*>>(arg.func.get());
	}
	else if (token == "asin")
	{
		result.func = std::make_shared<functions::ASinus<functions::Abstract*>>(arg.func.get());
	}
	else if (token == "acos")
	{
		result.func = std::make_shared<functions::ACosinus<functions::Abstract*>>(arg.func.get());
	}
	else if (token == "atg")
	{
		result.func = std::make_shared<functions::ATangent<functions::Abstract*>>(arg.func.get());
	}
	else if (token == "actg")
	{
		result.func = std::make_shared<functions::ACotangent<functions::Abstract*>>(arg.func.get());
	}
	else if (token == "sqrt")
	{
		result.func = std::make_shared<functions::Power<functions::Abstract*>>(arg.func.get(), 0.5);
	}
	else if (token == "sqr")
	{
		result.func = std::make_shared<functions::Power<functions::Abstract*>>(arg.func.get(), 2.0);
	}
	else if (token == "dx")
	{
		result.func = std::make_shared<Derivative<functions::Abstract*>>(arg.func.get());
	}
	else
	{
		throw std::runtime_error("Unknown unary function: " + token);
	}

	keep(result, arg);
	return result;
}

Expression Parser::makeBinary(const std::string& token, const Expression& left, const Expression& right)
{
	Expression result;
	result.token = token;

	if (token == "+")
	{
		result.func = std::make_shared<operations::Add<functions::Abstract*, functions::Abstract*>>(left.func.get(), right.func.get());
	}
	else if (token == "-")
	{
		result.func = std::make_shared<operations::Subtract<functions::Abstract*, functions::Abstract*>>(left.func.get(), right.func.get());
	}
	else if (token == "*")
	{
		result.func = std::make_shared<operations::Multiply<functions::Abstract*, functions::Abstract*>>(left.func.get(), right.func.get());
	}
	else if (token == "/")
	{
		result.func = std::make_shared<operations::Divide<functions::Abstract*, functions::Abstract*>>(left.func.get(), right.func.get());
	}
	else if (token == "^")
	{
		result.func = std::make_shared<functions::Exponent_Power<functions::Abstract*, functions::Abstract*>>(left.func.get(), right.func.get());
	}
	else
	{
		throw std::runtime_error("Unknown binary operator: " + token);
	}

	keep(result, left);
	keep(result, right);
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
			Expression result("pow", std::make_shared<functions::Power<functions::Abstract*>>(args[0].func.get(), args[1].number));
			keep(result, args[0]);
			return result;
		}
		return makeBinary("^", args[0], args[1]);
	}

	if (id == "piecewise")
	{
		if (args.size() < 3 || args.size() % 2 == 0)
		{
			throw std::runtime_error("piecewise expects: piecewise(border1, expr1, ..., lastExpr)");
		}

		std::vector<double> borders;
		std::vector<std::shared_ptr<functions::Abstract>> pieces;
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

		Expression result("piecewise", std::make_shared<functions::PiecewiseRuntime>(borders, pieces));
		for (const Expression& arg : args)
		{
			keep(result, arg);
		}
		return result;
	}

	if (id == "exp")
	{
		if (args.size() == 1)
		{
			Expression result("exp", std::make_shared<functions::Exponent<functions::Abstract*>>(std::exp(1.0), args[0].func.get()));
			keep(result, args[0]);
			return result;
		}
		if (args.size() == 2 && args[0].isNumber)
		{
			Expression result("exp", std::make_shared<functions::Exponent<functions::Abstract*>>(args[0].number, args[1].func.get()));
			keep(result, args[1]);
			return result;
		}
		throw std::runtime_error("exp expects exp(expr) or exp(number, expr)");
	}

	if (id == "log")
	{
		if (args.size() == 1)
		{
			auto base = std::make_shared<functions::Const>(std::exp(1.0));
			Expression result("log", std::make_shared<functions::Logarithm<functions::Abstract*, functions::Abstract*>>(base.get(), args[0].func.get()));
			keep(result, base);
			keep(result, args[0]);
			return result;
		}
		if (args.size() == 2)
		{
			Expression result("log", std::make_shared<functions::Logarithm<functions::Abstract*, functions::Abstract*>>(args[0].func.get(), args[1].func.get()));
			keep(result, args[0]);
			keep(result, args[1]);
			return result;
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
	return parsePrimary();
}

Expression Parser::parsePower()
{
	Expression left = parseUnary();
	if (consume('^'))
	{
		Expression right = parsePower();
		return makeBinary("^", left, right);
	}
	return left;
}

Expression Parser::parseMulDiv()
{
	Expression left = parsePower();
	for (;;)
	{
		if (consume('*'))
		{
			left = makeBinary("*", left, parsePower());
		}
		else if (consume('/'))
		{
			left = makeBinary("/", left, parsePower());
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

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double eval(const Expression& e) 
{
	return e.func->operator()(inputx);
}

double eval(const Expression& e, double x)
{
	return e.func->operator()(x);
}
