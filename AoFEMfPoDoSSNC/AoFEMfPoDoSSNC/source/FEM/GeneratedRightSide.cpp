#include "GeneratedRightSide.h"

#include "SymbolicExpression.h"

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
	const symbolic::SymbolicExpression pExpression = symbolic::parseExpression(p.c_str());
	const symbolic::SymbolicExpression qExpression = symbolic::parseExpression(q.c_str());
	const symbolic::SymbolicExpression uExpression = symbolic::parseExpression(exactU.c_str());

	const symbolic::SymbolicExpression uDerivative = symbolic::derivative(uExpression);
	const symbolic::SymbolicExpression flux = symbolic::binary("*", pExpression, uDerivative);
	const symbolic::SymbolicExpression fluxDerivative = symbolic::derivative(flux);
	const symbolic::SymbolicExpression reaction = symbolic::binary("*", qExpression, uExpression);

	return symbolic::toSimplifiedString(symbolic::binary("+", symbolic::unary("-", fluxDerivative), reaction));
}

