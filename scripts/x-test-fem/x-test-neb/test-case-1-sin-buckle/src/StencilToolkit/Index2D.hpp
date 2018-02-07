#pragma once

#include "Node2D.hpp"
#include "BlockRange.hpp"
#include "LocalBoundary2D.h"
#include "RemoteBoundary2D.h"
#include "SerialRange.hpp"
#include "Uncopyable.hpp"

namespace StencilToolkit
{

class Index2D : private Uncopyable
{
public:
	Index2D(index_t dim_x, index_t dim_y)
	: _node(nullptr)
	{
		_rx = new SerialRange(dim_x);
		_ry = new SerialRange(dim_y);
	}

	Index2D(Node2D &node, index_t dim_x, index_t dim_y)
	: _node(&node)
	{
		if (_node->size_x() == 1) _rx = new SerialRange(dim_x);
		else                      _rx = new BlockRange(dim_x, _node->size_x(), _node->rank_x());

		if (_node->size_y() == 1) _ry = new SerialRange(dim_y);
		else                      _ry = new BlockRange(dim_y, _node->size_y(), _node->rank_y());
	}

	virtual ~Index2D()
	{
		delete _rx;
		_rx = nullptr;
		delete _ry;
		_ry = nullptr;
	}

	virtual Range &get_range_x()
	{
		return *_rx;
	}

	virtual Range &get_range_y()
	{
		return *_ry;
	}

	Boundary *alloc_boundary(void *buffer, index_t type_bsize, index_t bound_x, index_t bound_y)
	{
		index_t dim_x = _rx->end() - _rx->begin();
		index_t dim_y = _ry->end() - _ry->begin();
		if (_node == nullptr)
		{
			return new LocalBoundary2D(buffer, type_bsize, bound_x, bound_y, dim_x, dim_y);
		}
		else
		{
			return new RemoteBoundary2D(_node, buffer, type_bsize, bound_x, bound_y, dim_x, dim_y);
		}
	}

private:
	// FIXME not good design
	Node2D *_node;

	Range *_rx;
	Range *_ry;

	Index2D();
};

}
