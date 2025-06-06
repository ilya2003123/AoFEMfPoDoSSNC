#pragma once

#include "../Functions/Abstract.h"

namespace operations  // �� ��� ��� �������� ������ namespace, ������������ ��������� ����� ������
{
	template<typename F1, typename F2>  // ��� ������� ��� 2 ������! �� ��� ��, �� ����������� ������� 2 ����
	class Add : public functions::Abstract // ����� �����
	{
	public:  // ������� ��� ���, ������� ��������� �� �������� ���� �����, ������� � ������� � ���� �����
		typedef Add<F1, F2> Type;
		Add(const F1& f1, const F2& f2) 
			:m_f1(f1), m_f2(f2)
		{
		}

		double operator()(double x) override  // ��, � ������ ����������� �� ��� �������
		{
			double f1 = 0;
			double f2 = 0;
			if constexpr (std::is_pointer_v<F1>)
				f1 += (*m_f1)(x);
			else
				f1 += m_f1(x);

			if constexpr (std::is_pointer_v<F2>)
				f2 += (*m_f2)(x);
			else
				f2 += m_f2(x);

			return f1 + f2;
		}

		F1 m_f1; // ���� ���������� ������ ����, ������ �������
		F2 m_f2;

	};
}