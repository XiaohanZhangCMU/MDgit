#pragma once

#include "Node1D.hpp"
#include "BlockRange.hpp"
#include "LocalBoundary1D.h"
#include "RemoteBoundary1D.h"
#include "SerialRange.hpp"
#include "Uncopyable.hpp"

namespace StencilToolkit
{

class Index1D : private Uncopyable
{
public:
	Index1D(index_t dim_size)
	: _node(nullptr)
	{
		_r = new SerialRange(dim_size);
	}

	Index1D(Node1D &node, index_t dim_size)
	: _node(&node)
	{
		if (_node->size() == 1) _r = new SerialRange(dim_size);
		else                    _r = new BlockRange(dim_size, _node->size(), _node->rank());
	}

	virtual ~Index1D()
	{
		delete _r;
		_r = nullptr;
	}

	virtual Range &get_range()
	{
		return *_r;
	}

	Boundary *alloc_boundary(void *buffer, index_t type_bsize, index_t bound)
	{
		index_t dim = _r->end() - _r->begin();
		if (_node == nullptr)
		{
			return new LocalBoundary1D(buffer, type_bsize, bound, dim);
		}
		else
		{
			return new RemoteBoundary1D(_node, buffer, type_bsize, bound, dim);
		}
	}

private:
	Node1D* _node;

	Range *_r;

	Index1D();
};

}
