#include "ConsoleInput.h"

#include <algorithm>
#include <cmath>
#include <cctype>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <vector>

namespace
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

	std::string toString(const SymbolicExpression& expression, int parentPrecedence = 0)
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

	struct SymbolicTerm
	{
		double coefficient = 1.0;
		std::map<std::string, int> powers;
		std::map<std::string, SymbolicExpression> bases;
	};

	bool nearlyZero(double value)
	{
		return std::abs(value) < 1e-10;
	}

	bool isInteger(double value)
	{
		return std::abs(value - std::round(value)) < 1e-10;
	}

	std::vector<std::string> sortedFactorKeys(const SymbolicTerm& term)
	{
		std::vector<std::string> keys;
		for (const auto& item : term.powers)
		{
			if (item.second != 0)
			{
				keys.push_back(item.first);
			}
		}

		std::sort(keys.begin(), keys.end(), [](const std::string& left, const std::string& right)
		{
			if (left == "x" && right != "x") return true;
			if (right == "x" && left != "x") return false;
			return left < right;
		});
		return keys;
	}

	void addFactor(SymbolicTerm& term, const SymbolicExpression& factor, int power = 1)
	{
		if (power == 0)
		{
			return;
		}
		if (factor.kind == SymbolicExpression::Kind::Number)
		{
			term.coefficient *= std::pow(factor.number, power);
			return;
		}
		if (factor.kind == SymbolicExpression::Kind::Unary && factor.token == "-")
		{
			if (power % 2 != 0)
			{
				term.coefficient *= -1.0;
			}
			addFactor(term, factor.args[0], power);
			return;
		}
		if (factor.kind == SymbolicExpression::Kind::Binary && factor.token == "*" && power == 1)
		{
			addFactor(term, factor.args[0]);
			addFactor(term, factor.args[1]);
			return;
		}
		if (factor.kind == SymbolicExpression::Kind::Binary
			&& factor.token == "^"
			&& factor.args[1].kind == SymbolicExpression::Kind::Number
			&& isInteger(factor.args[1].number))
		{
			addFactor(term, factor.args[0], power * static_cast<int>(std::round(factor.args[1].number)));
			return;
		}

		const std::string key = toString(factor);
		term.powers[key] += power;
		term.bases.emplace(key, factor);
	}

	SymbolicTerm multiplyTerms(SymbolicTerm left, const SymbolicTerm& right)
	{
		left.coefficient *= right.coefficient;
		for (const auto& item : right.powers)
		{
			left.powers[item.first] += item.second;
		}
		for (const auto& item : right.bases)
		{
			left.bases.emplace(item.first, item.second);
		}
		return left;
	}

	std::vector<SymbolicTerm> collectTerms(const SymbolicExpression& expression)
	{
		switch (expression.kind)
		{
		case SymbolicExpression::Kind::Number:
			return { SymbolicTerm{ expression.number, {}, {} } };
		case SymbolicExpression::Kind::Variable:
		case SymbolicExpression::Kind::Piecewise:
		{
			SymbolicTerm term;
			addFactor(term, expression);
			return { term };
		}
		case SymbolicExpression::Kind::Unary:
			if (expression.token == "-")
			{
				std::vector<SymbolicTerm> terms = collectTerms(expression.args[0]);
				for (SymbolicTerm& term : terms)
				{
					term.coefficient *= -1.0;
				}
				return terms;
			}
			{
				SymbolicTerm term;
				addFactor(term, expression);
				return { term };
			}
		case SymbolicExpression::Kind::Binary:
			if (expression.token == "+")
			{
				std::vector<SymbolicTerm> terms = collectTerms(expression.args[0]);
				std::vector<SymbolicTerm> right = collectTerms(expression.args[1]);
				terms.insert(terms.end(), right.begin(), right.end());
				return terms;
			}
			if (expression.token == "-")
			{
				std::vector<SymbolicTerm> terms = collectTerms(expression.args[0]);
				std::vector<SymbolicTerm> right = collectTerms(expression.args[1]);
				for (SymbolicTerm& term : right)
				{
					term.coefficient *= -1.0;
				}
				terms.insert(terms.end(), right.begin(), right.end());
				return terms;
			}
			if (expression.token == "*")
			{
				std::vector<SymbolicTerm> result;
				const std::vector<SymbolicTerm> left = collectTerms(expression.args[0]);
				const std::vector<SymbolicTerm> right = collectTerms(expression.args[1]);
				for (const SymbolicTerm& leftTerm : left)
				{
					for (const SymbolicTerm& rightTerm : right)
					{
						result.push_back(multiplyTerms(leftTerm, rightTerm));
					}
				}
				return result;
			}
			if (expression.token == "/" && expression.args[1].kind == SymbolicExpression::Kind::Number)
			{
				std::vector<SymbolicTerm> terms = collectTerms(expression.args[0]);
				for (SymbolicTerm& term : terms)
				{
					term.coefficient /= expression.args[1].number;
				}
				return terms;
			}
			{
				SymbolicTerm term;
				addFactor(term, expression);
				return { term };
			}
		}

		return {};
	}

	std::string termKey(const SymbolicTerm& term)
	{
		std::string key;
		for (const std::string& factor : sortedFactorKeys(term))
		{
			key += factor + "^" + std::to_string(term.powers.at(factor)) + "|";
		}
		return key;
	}

	int xDegree(const SymbolicTerm& term)
	{
		const auto found = term.powers.find("x");
		return found == term.powers.end() ? 0 : found->second;
	}

	std::string factorToString(const std::string& key, int power)
	{
		if (power == 1)
		{
			return key;
		}
		if (key.find_first_of("+-*/") != std::string::npos)
		{
			return "(" + key + ")^" + std::to_string(power);
		}
		return key + "^" + std::to_string(power);
	}

	std::string termBodyToString(const SymbolicTerm& term)
	{
		const double absCoefficient = std::abs(term.coefficient);
		const std::vector<std::string> factors = sortedFactorKeys(term);

		if (factors.empty())
		{
			return formatNumber(absCoefficient);
		}

		std::string result;
		if (std::abs(absCoefficient - 1.0) > 1e-10)
		{
			result = formatNumber(absCoefficient);
		}

		for (const std::string& factor : factors)
		{
			if (!result.empty())
			{
				result += "*";
			}
			result += factorToString(factor, term.powers.at(factor));
		}
		return result;
	}

	std::string toSimplifiedString(const SymbolicExpression& expression)
	{
		std::map<std::string, SymbolicTerm> combined;
		for (const SymbolicTerm& term : collectTerms(expression))
		{
			if (nearlyZero(term.coefficient))
			{
				continue;
			}

			const std::string key = termKey(term);
			auto found = combined.find(key);
			if (found == combined.end())
			{
				combined.emplace(key, term);
			}
			else
			{
				found->second.coefficient += term.coefficient;
			}
		}

		std::vector<SymbolicTerm> terms;
		for (const auto& item : combined)
		{
			if (!nearlyZero(item.second.coefficient))
			{
				terms.push_back(item.second);
			}
		}

		if (terms.empty())
		{
			return "0";
		}

		std::sort(terms.begin(), terms.end(), [](const SymbolicTerm& left, const SymbolicTerm& right)
		{
			const std::vector<std::string> leftKeys = sortedFactorKeys(left);
			const std::vector<std::string> rightKeys = sortedFactorKeys(right);
			const bool leftConstant = leftKeys.empty();
			const bool rightConstant = rightKeys.empty();
			if (leftConstant != rightConstant) return leftConstant;

			const bool leftPurePower = leftKeys.size() == 1 && leftKeys[0] == "x";
			const bool rightPurePower = rightKeys.size() == 1 && rightKeys[0] == "x";
			if (leftPurePower != rightPurePower) return leftPurePower;

			const int leftDegree = xDegree(left);
			const int rightDegree = xDegree(right);
			if (leftPurePower && rightPurePower && leftDegree != rightDegree) return leftDegree < rightDegree;

			return termKey(left) < termKey(right);
		});

		std::string result;
		for (const SymbolicTerm& term : terms)
		{
			const bool negative = term.coefficient < 0.0;
			const std::string body = termBodyToString(term);

			if (result.empty())
			{
				result = negative ? "-" + body : body;
			}
			else
			{
				result += negative ? "-" : "+";
				result += body;
			}
		}
		return result;
	}

	class SymbolicParser
	{
	public:
		explicit SymbolicParser(const char* input)
			: m_input(input)
		{
		}

		SymbolicExpression parse()
		{
			SymbolicExpression result = parseAddSub();
			skipSpaces();
			if (*m_input != '\0')
			{
				throw std::runtime_error("Unexpected tail in generated f display expression");
			}
			return result;
		}

	private:
		void skipSpaces()
		{
			while (std::isspace(static_cast<unsigned char>(*m_input)))
			{
				++m_input;
			}
		}

		bool consume(char expected)
		{
			skipSpaces();
			if (*m_input != expected)
			{
				return false;
			}
			++m_input;
			return true;
		}

		std::string parseIdentifier()
		{
			skipSpaces();
			std::string result;
			while (std::isalpha(static_cast<unsigned char>(*m_input)) || *m_input == '_')
			{
				result.push_back(*m_input++);
			}
			return result;
		}

		SymbolicExpression parseNumber()
		{
			skipSpaces();
			char* end = nullptr;
			const double value = std::strtod(m_input, &end);
			if (end == m_input)
			{
				throw std::runtime_error("Expected number");
			}
			m_input = end;
			return number(value);
		}

		std::vector<SymbolicExpression> parseArguments()
		{
			std::vector<SymbolicExpression> args;
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
					throw std::runtime_error("Expected ',' or ')' in function arguments");
				}
			}
		}

		SymbolicExpression parsePrimary()
		{
			skipSpaces();
			if (std::isdigit(static_cast<unsigned char>(*m_input)) || *m_input == '.')
			{
				return parseNumber();
			}
			if (consume('('))
			{
				SymbolicExpression inside = parseAddSub();
				if (!consume(')'))
				{
					throw std::runtime_error("Expected ')'");
				}
				return inside;
			}

			const std::string id = parseIdentifier();
			if (id.empty())
			{
				throw std::runtime_error("Unexpected character in expression");
			}
			if (id == "x")
			{
				return variable();
			}
			if (id == "pi")
			{
				return number(std::acos(-1.0));
			}
			if (id == "e")
			{
				return number(std::exp(1.0));
			}

			std::vector<SymbolicExpression> args;
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
				return binary("^", args[0], args[1]);
			}
			if (id == "pw" || id == "piecewise")
			{
				if (args.size() < 3 || args.size() % 2 == 0)
				{
					throw std::runtime_error("pw expects: pw(border1, expr1, ..., lastExpr)");
				}
				SymbolicExpression expression;
				expression.kind = SymbolicExpression::Kind::Piecewise;
				expression.token = "pw";
				expression.args = std::move(args);
				return expression;
			}
			if (id == "exp" && args.size() == 2)
			{
				return binary("^", args[0], args[1]);
			}
			if (args.size() != 1)
			{
				throw std::runtime_error(id + " expects 1 argument");
			}
			return unary(id, args[0]);
		}

		SymbolicExpression parsePower()
		{
			SymbolicExpression left = parsePrimary();
			if (consume('^'))
			{
				return binary("^", left, parseUnary());
			}
			return left;
		}

		SymbolicExpression parseUnary()
		{
			if (consume('+'))
			{
				return parseUnary();
			}
			if (consume('-'))
			{
				return unary("-", parseUnary());
			}
			return parsePower();
		}

		SymbolicExpression parseMulDiv()
		{
			SymbolicExpression left = parseUnary();
			for (;;)
			{
				if (consume('*'))
				{
					left = binary("*", left, parseUnary());
				}
				else if (consume('/'))
				{
					left = binary("/", left, parseUnary());
				}
				else
				{
					return left;
				}
			}
		}

		SymbolicExpression parseAddSub()
		{
			SymbolicExpression left = parseMulDiv();
			for (;;)
			{
				if (consume('+'))
				{
					left = binary("+", left, parseMulDiv());
				}
				else if (consume('-'))
				{
					left = binary("-", left, parseMulDiv());
				}
				else
				{
					return left;
				}
			}
		}

		const char* m_input;
	};
}

void printExpressionHelp()
{
	std::cout << "Expression syntax: use x, +, -, *, /, ^ and functions like sin(...), cos(...), exp(...), log(...).\n";
	std::cout << "Piecewise syntax: pw(border1, expr1, border2, expr2, lastExpr).\n";
	std::cout << "The longer name piecewise(...) is also supported.\n";
	std::cout << "It means expr1 for x <= border1, expr2 for x <= border2, otherwise lastExpr.\n\n";
}

ProblemDefinition problemFromConsole()
{
	ProblemDefinition problem;
	printExpressionHelp();

	std::cout << "N: ";
	std::cin >> problem.intervals;

	std::cout << "l: ";
	std::cin >> problem.length;

	std::cout << "m: ";
	std::cin >> problem.limiter;
	std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

	std::cout << "p(x): ";
	std::getline(std::cin, problem.p);

	std::cout << "q(x): ";
	std::getline(std::cin, problem.q);

	std::cout << "f(x): ";
	std::getline(std::cin, problem.f);

	std::cout << "Exact u(x), empty if unknown: ";
	std::getline(std::cin, problem.exactU);

	return problem;
}

std::string buildGeneratedRightSideExpression(
	const std::string& p,
	const std::string& q,
	const std::string& exactU)
{
	return "-(dx((" + p + ")*dx(" + exactU + ")))+(" + q + ")*(" + exactU + ")";
}

std::string buildGeneratedRightSideDisplayExpression(
	const std::string& p,
	const std::string& q,
	const std::string& exactU)
{
	const SymbolicExpression pExpression = SymbolicParser(p.c_str()).parse();
	const SymbolicExpression qExpression = SymbolicParser(q.c_str()).parse();
	const SymbolicExpression uExpression = SymbolicParser(exactU.c_str()).parse();

	const SymbolicExpression uDerivative = derivative(uExpression);
	const SymbolicExpression flux = binary("*", pExpression, uDerivative);
	const SymbolicExpression fluxDerivative = derivative(flux);
	const SymbolicExpression reaction = binary("*", qExpression, uExpression);

	return toSimplifiedString(binary("+", unary("-", fluxDerivative), reaction));
}

ProblemDefinition problemFromExactSolution()
{
	ProblemDefinition problem;
	problem.generatedRightSide = true;
	problem.name = "Generated from exact solution";
	printExpressionHelp();

	std::cout << "N: ";
	std::cin >> problem.intervals;

	std::cout << "l: ";
	std::cin >> problem.length;

	std::cout << "m: ";
	std::cin >> problem.limiter;
	std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

	std::cout << "p(x): ";
	std::getline(std::cin, problem.p);

	std::cout << "q(x): ";
	std::getline(std::cin, problem.q);

	std::cout << "Exact u(x): ";
	std::getline(std::cin, problem.exactU);

	problem.f = buildGeneratedRightSideExpression(problem.p, problem.q, problem.exactU);
	return problem;
}

int chooseRunMode()
{
	std::cout << "Choose run mode:\n";
	std::cout << "  1 - enter p(x), q(x), f(x) manually\n";
	std::cout << "  2 - enter exact u(x), p(x), q(x); generate f(x)\n";
	std::cout << "Mode: ";

	int mode = 1;
	std::cin >> mode;
	std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

	if (mode != 1 && mode != 2)
	{
		throw std::runtime_error("Unknown run mode. Use 1 or 2.");
	}

	return mode;
}
