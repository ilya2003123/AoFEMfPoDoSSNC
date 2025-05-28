#pragma once

#include <assert.h>
#include "../Derivative/derivative.h"
#include "../Functions/functions.h"
#include "../Operations/operations.h"

using namespace utils;

void test()
{
	const double base = 2.72;
	auto pow1 = (Pow(X, 2) / 4);
	auto exp = Exp(base, pow1);
	auto cos = Cos(X);
	auto dcos = derivative(cos);

	std::string str1 = "exp(" + std::to_string(base) + ", (pow(x, 2) / 4))";
	inputx = 1;

	Parser p1(str1.c_str());
	auto q1 = p1.parse();
	auto result1 = eval(q1);


	auto DIV = functions::Const(3) / X;
	auto dDIV = derivative(DIV);
	std::cout << dDIV(3) << std::endl;
	// std::cout << result1 << std::endl;

	assert(result1 - exp(inputx) <= 1e-4);

	assert(Sin(0) == 0);
	assert(Cos(0) == 1);
	assert(Pow(2, 2) == 4);
	assert(Exp(2, 4) == 16);
	assert(derivative(Exp(2,X))(1) - 1.38629 <= 1e-5);
	assert(abs(exp(1) - 1.2842) <= 1e-4);
	assert(dcos(10000) - 0.30561 <= 1e-5);
	assert(derivative(Log(X, X))(2) == 0);
	assert(derivative(ATg(X))(1) == 0.5);
	assert(derivative(ACos(X))(0.5) + 1.15470 <= 1e-5);
	assert(derivative(ASin(X))(0.5) - 1.15470 <= 1e-5);
	assert(derivative(ACtg(X))(1) == -0.5);
	assert(derivative(X + Cos(X))(0) == 1);
	assert(derivative(Pow(X, 3))(2) == 12);
	assert(derivative(X * Cos(X))(3) + 1.41335 <= 1e-5);

	std::cout << "All test done" << std::endl;
}

void exmain()
{
	std::string str2 = "dx(cos(x))";
	Parser p2(str2.c_str());

	inputx = 10000;

	auto q2 = p2.parse();
	auto result2 = eval(q2);
	std::cout << result2 << std::endl;

	//std::cout << dcos(1);

	try
	{
		std::string str;
		std::cout << "Function: ";
		std::getline(std::cin, str);
		std::cout << std::endl;
		std::cout << "X: ";
		std::cin >> inputx;
		std::cout << std::endl;

		Parser p(str.c_str());
		auto q = p.parse();
		auto result = eval(q);
		std::cout << result << std::endl;
	}
	catch (std::runtime_error& Error) { std::cout << Error.what() << std::endl; }
}