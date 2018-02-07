#pragma once

#include "Node.hpp"

namespace StencilToolkit
{

class Node3D : public Node
{
public:
	Node3D(index_t x, index_t y, index_t z)
	{
		_size_x = x;
		_size_y = y;
		_size_z = z;

		index_t size = MachineInfo::get_instance().get_size();
		if (size != (_size_x * _size_y * _size_z))
		{
			MachineInfo::error("error on Node3D constructor: wrong size");
		}

		index_t rank = MachineInfo::get_instance().get_rank();

		_rank_x = rank % _size_x;
		_rank_y = (rank / _size_x) % _size_y;
		_rank_z = rank / (_size_x * _size_y);
	}

	index_t operator()(index_t x, index_t y, index_t z)
	{
		x = x % _size_x;
		index_t temp_x = (_rank_x + x + _size_x) % _size_x;
		y = y % _size_y;
		index_t temp_y = (_rank_y + y + _size_y) % _size_y;
		z = z % _size_z;
		index_t temp_z = (_rank_z + z + _size_z) % _size_z;
		return calc_linear_rank(temp_x, temp_y, temp_z);
	}

	index_t size_x()
	{
		return _size_x;
	}

	index_t size_y()
	{
		return _size_y;
	}

	index_t size_z()
	{
		return _size_z;
	}

	index_t rank_x()
	{
		return _rank_x;
	}

	index_t rank_y()
	{
		return _rank_y;
	}

	index_t rank_z()
	{
		return _rank_z;
	}

	virtual bool is_master()
	{
		return (_rank_x == 0) && (_rank_y == 0) && (_rank_z == 0);
	}

private:
	index_t calc_linear_rank(index_t x, index_t y, index_t z)
	{
		return x + (y * _size_x) + (z * _size_x * _size_y);
	}

	index_t _size_x;
	index_t _size_y;
	index_t _size_z;

	index_t _rank_x;
	index_t _rank_y;
	index_t _rank_z;

	Node3D();
};

}
