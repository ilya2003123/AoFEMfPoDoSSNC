#include "SymbolicExpression.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <string>
#include <vector>

namespace symbolic
{
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
}
