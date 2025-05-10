#pragma warning(disable: 26451)

#include <iostream>
#include <string>
#include <vector>
#include "Functions/functions.h"
#include "Operations/operations.h"
#include "Utils/Utils.h"
#include "Derivative/derivative.h"
#include "Parser/Parser.h"
#include "test/test.h"
//#include <iomanip> - для вывода, можно отрегулировать сколько знаков после запятой выводится!

using namespace utils;  // Просто потому что надо! :D


int main()
{
	test();

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

	return 0;
}