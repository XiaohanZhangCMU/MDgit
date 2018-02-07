#pragma once

#include "Node.hpp"

namespace StencilToolkit
{

class Node1D : public Node
{
public:
	Node1D(index_t x)
	{
		_size = x;

		index_t size = MachineInfo::get_instance().get_size();
		if (size != _size)
		{
			MachineInfo::error("error on Node1D constructor: wrong size");
		}

		_rank = MachineInfo::get_instance().get_rank();
	}

	index_t operator()(index_t x)
	{
		x = x % _size;
		return (_rank + x + _size) % _size;
	}

	index_t size()
	{
		return _size;
	}

	index_t rank()
	{
		return _rank;
	}

	virtual bool is_master()
	{
		return (_rank == 0);
	}

private:
	index_t _size;
	index_t _rank;

	Node1D();
};

}
