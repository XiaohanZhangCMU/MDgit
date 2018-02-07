#pragma once

#include "common.h"
#include "Boundary.hpp"
#include "Node3D.hpp"
#include "Range.hpp"

namespace StencilToolkit
{

class RemoteBoundary3D : public Boundary
{
public:
	RemoteBoundary3D(Node3D *node, void *buffer, index_t type_bsize,
	                 index_t bound_x, index_t bound_y, index_t bound_z,
	                 index_t dim_x, index_t dim_y, index_t dim_z)
	: _node(node), _buffer(buffer), _type_bsize(type_bsize),
	  _boundary_size_x(bound_x), _boundary_size_y(bound_y), _boundary_size_z(bound_z),
	  _dim_size_x(dim_x), _dim_size_y(dim_y), _dim_size_z(dim_z)
	{
		init_type();
	}

	~RemoteBoundary3D()
	{
		finalize_type();
	}

	virtual void do_sync();

private:
	index_t calc_bpos(index_t type_bsize, index_t x, index_t y, index_t z,
	                  index_t start_index, index_t diff_y, index_t diff_z)
	{
		return type_bsize * (start_index + x + (y * diff_y) + (z * diff_z));
	}

	void init_type();
	void finalize_type();

	Node3D *_node;

	void *_buffer;
	index_t _type_bsize;

	index_t _boundary_size_x;
	index_t _boundary_size_y;
	index_t _boundary_size_z;

	index_t _dim_size_x;
	index_t _dim_size_y;
	index_t _dim_size_z;

	void *_boundary_type;

	RemoteBoundary3D();
};

}
