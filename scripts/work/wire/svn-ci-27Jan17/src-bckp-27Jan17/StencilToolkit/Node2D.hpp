#pragma once

#include "Node.hpp"

namespace StencilToolkit
{

class Node2D : public Node
{
public:
	Node2D(index_t x, index_t y)
	{
		_size_x = x;
		_size_y = y;

		index_t size = MachineInfo::get_instance().get_size();
		if (size != (_size_x * _size_y))
		{
			MachineInfo::error("error on Node2D constructor: wrong size");
		}

		index_t rank = MachineInfo::get_instance().get_rank();

		_rank_x = rank % _size_x;
		_rank_y = rank / _size_x;
	}

	index_t operator()(index_t x, index_t y)
	{
		x = x % _size_x;
		index_t temp_x = (_rank_x + x + _size_x) % _size_x;
		y = y % _size_y;
		index_t temp_y = (_rank_y + y + _size_y) % _size_y;
		return calc_linear_rank(temp_x, temp_y);
	}

	index_t size_x()
	{
		return _size_x;
	}

	index_t size_y()
	{
		return _size_y;
	}

	index_t rank_x()
	{
		return _rank_x;
	}

	index_t rank_y()
	{
		return _rank_y;
	}

	virtual bool is_master()
	{
		return (_rank_x == 0) && (_rank_y == 0);
	}

private:
	index_t calc_linear_rank(index_t x, index_t y)
	{
		return x + (y * _size_x);
	}

	index_t _size_x;
	index_t _size_y;

	index_t _rank_x;
	index_t _rank_y;

	Node2D();
};

}
