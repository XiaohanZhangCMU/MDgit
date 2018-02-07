#pragma once

#include "common.h"
#include "Boundary.hpp"
#include "Node2D.hpp"
#include "Range.hpp"

namespace StencilToolkit
{

class RemoteBoundary2D : public Boundary
{
public:
	RemoteBoundary2D(Node2D *node, void *buffer, index_t type_bsize,
	                 index_t bound_x, index_t bound_y,
	                 index_t dim_x, index_t dim_y)
	: _node(node), _buffer(buffer), _type_bsize(type_bsize),
	  _boundary_size_x(bound_x), _boundary_size_y(bound_y),
	  _dim_size_x(dim_x), _dim_size_y(dim_y)
	{
		init_type();
	}

	~RemoteBoundary2D()
	{
		finalize_type();
	}

	virtual void do_sync();

private:
	index_t calc_bpos(index_t type_bsize, index_t x, index_t y, index_t start_index, index_t diff_y)
	{
		return type_bsize * (start_index + x + (y * diff_y));
	}

	void init_type();
	void finalize_type();

	Node2D *_node;

	void *_buffer;
	index_t _type_bsize;

	index_t _boundary_size_x;
	index_t _boundary_size_y;

	index_t _dim_size_x;
	index_t _dim_size_y;

	void *_boundary_type;

	RemoteBoundary2D();
};

}
