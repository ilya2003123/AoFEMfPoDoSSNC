#pragma once

#include "../Functions/functions.h"
#include <memory>
#include <string>
#include <vector>

struct Expression
{
	Expression() = default;
	Expression(std::string token, std::shared_ptr<functions::Abstract> f);

	std::string token;
	std::shared_ptr<functions::Abstract> func;
	double number = 0.0;
	bool isNumber = false;
};

class Parser
{
public:
	explicit Parser(const char* input);
	Expression parse();

private:
	void skipSpaces();
	bool match(char expected);
	bool consume(char expected);
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

	const char* input;
	const char* begin;
};

double eval(const Expression& e, double x);
