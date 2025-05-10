#pragma once

#pragma warning(disable: 26495)

#include "../Functions/functions.h"
#include "../Operations/operations.h"
#include "../Derivative/derivative.h"
#include <string>
#include <vector>
#include <cctype>
#include <cstring>
#include <stdexcept>

double inputx = 0;

struct Expression 
{
	Expression(std::string token, functions::Abstract* f = nullptr) : token(token), func(f) {}
	Expression(std::string token, Expression a) : token(token), args{ a } {}
	Expression(std::string token, Expression a, Expression b) : token(token), args{ a, b } {}
	Expression(std::string token, double base, functions::Abstract* f)
		: token(token), number(base), func(f) {}
	Expression(std::string token, functions::Abstract* f1, functions::Abstract* f2)
		: token(token), func(f1), func2(f2) {}

	std::string token;
	functions::Abstract* func;
	functions::Abstract* func2;
	std::vector<Expression> args;
	double number = 0;
};

class Parser 
{
public:
	explicit Parser(const char* input)
		: input(input) 
	{
	}
	Expression parse();

private:
	std::string parse_token();
	Expression parse_simple_expression();
	Expression parse_binary_expression(int min_priority);

	const char* input;
};

std::string Parser::parse_token() 
{
	while (std::isspace(*input)) ++input;

	if (std::isdigit(*input)) 
	{
		std::string number;
		while (std::isdigit(*input) || *input == '.') number.push_back(*input++);
		return number;
	}

	static const std::string tokens[] =
	{ "+", "-", "*", "/", "acos", "sin", "cos", "(", ")", "x", "actg", "asin", "atg", "ctg", "sqrt",
	"pow", "log", "exp", "sqr", "tg", ",", "dx", "^" };
	for (auto& t : tokens) 
	{
		if (std::strncmp(input, t.c_str(), t.size()) == 0) 
		{
			input += t.size();
			return t;
		}
	}

	return "";
}

Expression Parser::parse_simple_expression()
{
	auto token = parse_token();
	if (token.empty()) throw std::runtime_error("Invalid input");

	if (token == "(") {
		auto result = parse();
		if (parse_token() != ")") throw std::runtime_error("Expected ')'");
		return result;
	}

	if (std::isdigit(token[0]))
		return Expression(token, new functions::Const(stoi(token)));

	if (token == "x")
		return Expression(token, new functions::Simple);

	auto rightside = parse_simple_expression();

	if (token == "-")
		return Expression(token, new operations::UnarMinus(rightside.func));

	if (token == "sin")
		return Expression(token, new functions::Sinus(rightside.func));

	if (token == "cos")
		return Expression(token, new functions::Cosinus(rightside.func));

	if (token == "actg")
		return Expression(token, new functions::ACotangent(rightside.func));

	if (token == "acos")
		return Expression(token, new functions::ACosinus(rightside.func));

	if (token == "asin")
		return Expression(token, new functions::ASinus(rightside.func));

	if (token == "atg")
		return Expression(token, new functions::ATangent(rightside.func));

	if (token == "ctg")
		return Expression(token, new functions::Cotangent(rightside.func));

	if (token == "tg")
		return Expression(token, new functions::Tangent(rightside.func));

	if (token == "exp")
		return Expression(token, new functions::Exponent(rightside.number, rightside.func));

	if (token == "log")
		return Expression(token, new functions::Logarithm(rightside.func, rightside.func2));

	if (token == "pow")
		return Expression(token, new functions::Power(rightside.func, rightside.number));

	if (token == "sqrt")
		return Expression(token, new functions::Power(rightside.func, 0.5));

	if (token == "sqr")
		return Expression(token, new functions::Power(rightside.func, 2));

	if (token == "dx")
		return Expression(token, new Derivative(rightside.func));

	return Expression(token, rightside);
}

int get_priority(const std::string& binary_op) 
{
	if (binary_op == "+") return 1;
	if (binary_op == "-") return 1;
	if (binary_op == "*") return 2;
	if (binary_op == "/") return 2;
	if (binary_op == ",") return 3;
	if (binary_op == "^") return 4;
						  return 0;
}

Expression Parser::parse_binary_expression(int min_priority)
{
	auto left_expr = parse_simple_expression();

	for (;;) 
	{
		auto op = parse_token();
		auto priority = get_priority(op);
		if (priority <= min_priority) 
		{
			input -= op.size();
			return left_expr;
		}

		auto right_expr = parse_binary_expression(priority);

		if (op == "+")
			left_expr = Expression(op, new operations::Add(left_expr.func, right_expr.func));
		else if (op == "-")
			left_expr = Expression(op, new operations::Subtract(left_expr.func, right_expr.func));
		else if (op == "*")
			left_expr = Expression(op, new operations::Multiply(left_expr.func, right_expr.func));
		else if (op == "/")
			left_expr = Expression(op, new operations::Divide(left_expr.func, right_expr.func));
		else if (op == ",")
			left_expr = Expression(op, left_expr.func, right_expr.func);																
		else if (op == "^")
			left_expr = Expression(op, new functions::Exponent_Power(left_expr.func, right_expr.func));
		else
			left_expr = Expression(op, left_expr, right_expr);
	}
}

Expression Parser::parse() 
{
	return parse_binary_expression(0);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double eval(const Expression& e) 
{
	auto y = e.func->operator()(inputx);

	return y;
}