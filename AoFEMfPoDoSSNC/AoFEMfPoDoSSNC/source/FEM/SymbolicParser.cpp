#include "SymbolicExpression.h"

#include <cctype>
#include <cmath>
#include <cstdlib>
#include <stdexcept>
#include <string>
#include <vector>

namespace symbolic
{
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

	SymbolicExpression parseExpression(const char* input)
	{
		return SymbolicParser(input).parse();
	}
}

