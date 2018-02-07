#pragma once

#include "Node3D.hpp"
#include "BlockRange.hpp"
#include "LocalBoundary3D.h"
#include "RemoteBoundary3D.h"
#include "SerialRange.hpp"
#include "Uncopyable.hpp"

namespace StencilToolkit
{

class Index3D : private Uncopyable
{
public:
	Index3D(index_t dim_x, index_t dim_y, index_t dim_z)
	: _node(nullptr)
	{
		_rx = new SerialRange(dim_x);
		_ry = new SerialRange(dim_y);
		_rz = new SerialRange(dim_z);
	}

	Index3D(Node3D &node, index_t dim_x, index_t dim_y, index_t dim_z)
	: _node(&node)
	{
		if (_node->size_x() == 1) _rx = new SerialRange(dim_x);
		else                      _rx = new BlockRange(dim_x, _node->size_x(), _node->rank_x());

		if (_node->size_y() == 1) _ry = new SerialRange(dim_y);
		else                      _ry = new BlockRange(dim_y, _node->size_y(), _node->rank_y());

		if (_node->size_z() == 1) _rz = new SerialRange(dim_z);
		else                      _rz = new BlockRange(dim_z, _node->size_z(), _node->rank_z());
	}

	virtual ~Index3D()
	{
		delete _rx;
		_rx = nullptr;
		delete _ry;
		_ry = nullptr;
		delete _rz;
		_rz = nullptr;
	}

	virtual Range &get_range_x()
	{
		return *_rx;
	}

	virtual Range &get_range_y()
	{
		return *_ry;
	}

	virtual Range &get_range_z()
	{
		return *_rz;
	}

	Boundary *alloc_boundary(void *buffer, index_t type_bsize, index_t bound_x, index_t bound_y, index_t bound_z)
	{
		index_t dim_x = _rx->end() - _rx->begin();
		index_t dim_y = _ry->end() - _ry->begin();
		index_t dim_z = _rz->end() - _rz->begin();
		if (_node == nullptr)
		{
			return new LocalBoundary3D(buffer, type_bsize, bound_x, bound_y, bound_z, dim_x, dim_y, dim_z);
		}
		else
		{
			return new RemoteBoundary3D(_node, buffer, type_bsize, bound_x, bound_y, bound_z, dim_x, dim_y, dim_z);
		}
	}

private:
	Node3D *_node;

	Range *_rx;
	Range *_ry;
	Range *_rz;

	Index3D();
};

}
