#pragma once

#include <assert.h>
#include "../Derivative/derivative.h"
#include "../Functions/functions.h"
#include "../Operations/operations.h"
#include "../Parser/Parser.h"

void test()
{
	const double base = 2.72;
	auto pow1 = (utils::Pow(utils::X, 2) / 4);
	auto exp = utils::Exp(base, pow1);
	auto cos = utils::Cos(utils::X);
	auto dcos = derivative(cos);

	std::string str1 = "exp(" + std::to_string(base) + ", (pow(x, 2) / 4))";
	double x = 1;

	Parser p1(str1.c_str());
	auto q1 = p1.parse();
	auto result1 = eval(q1, x);


	auto DIV = functions::Const(3) / utils::X;
	auto dDIV = derivative(DIV);
	std::cout << dDIV(3) << std::endl;
	// std::cout << result1 << std::endl;

	assert(result1 - exp(x) <= 1e-4);

	assert(utils::Sin(0) == 0);
	assert(utils::Cos(0) == 1);
	assert(utils::Pow(2, 2) == 4);
	assert(utils::Exp(2, 4) == 16);
	assert(derivative(utils::Exp(2, utils::X))(1) - 1.38629 <= 1e-5);
	assert(abs(exp(1) - 1.2842) <= 1e-4);
	assert(dcos(10000) - 0.30561 <= 1e-5);
	assert(derivative(utils::Log(utils::X, utils::X))(2) == 0);
	assert(derivative(utils::ATg(utils::X))(1) == 0.5);
	assert(derivative(utils::ACos(utils::X))(0.5) + 1.15470 <= 1e-5);
	assert(derivative(utils::ASin(utils::X))(0.5) - 1.15470 <= 1e-5);
	assert(derivative(utils::ACtg(utils::X))(1) == -0.5);
	assert(derivative(utils::X + utils::Cos(utils::X))(0) == 1);
	assert(derivative(utils::Pow(utils::X, 3))(2) == 12);
	assert(derivative(utils::X * utils::Cos(utils::X))(3) + 1.41335 <= 1e-5);

	std::cout << "All test done" << std::endl;
}

void exmain()
{
	std::string str2 = "dx(cos(x))";
	Parser p2(str2.c_str());

	double x = 10000;

	auto q2 = p2.parse();
	auto result2 = eval(q2, x);
	std::cout << result2 << std::endl;

	//std::cout << dcos(1);

	try
	{
		std::string str;
		std::cout << "Function: ";
		std::getline(std::cin, str);
		std::cout << std::endl;
		std::cout << "X: ";
		std::cin >> x;
		std::cout << std::endl;

		Parser p(str.c_str());
		auto q = p.parse();
		auto result = eval(q, x);
		std::cout << result << std::endl;
	}
	catch (std::runtime_error& Error) { std::cout << Error.what() << std::endl; }
}
