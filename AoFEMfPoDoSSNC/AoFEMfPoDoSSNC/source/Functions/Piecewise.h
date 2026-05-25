#pragma once

#include "Abstract.h"
#include <memory>
#include <utility>
#include <vector>

namespace functions
{
	class Piecewise : public Abstract
	{
	public:
		Piecewise(std::vector<double> borders, std::vector<std::shared_ptr<Abstract>> pieces)
			: m_borders(std::move(borders)), m_pieces(std::move(pieces))
		{
		}

		double operator()(double x) override
		{
			for (int i = 0; i < static_cast<int>(m_borders.size()); ++i)
			{
				if (x <= m_borders[i] + 1e-12)
				{
					return (*m_pieces[i])(x);
				}
			}
			return (*m_pieces.back())(x);
		}

	private:
		std::vector<double> m_borders;
		std::vector<std::shared_ptr<Abstract>> m_pieces;
	};
}
