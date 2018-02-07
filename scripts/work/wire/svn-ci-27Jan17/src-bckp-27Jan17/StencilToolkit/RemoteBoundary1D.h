#pragma once

#include "common.h"
#include "Boundary.hpp"
#include "Node1D.hpp"
#include "Range.hpp"

namespace StencilToolkit
{

class RemoteBoundary1D : public Boundary
{
public:
	RemoteBoundary1D(Node1D *node, void *buffer, index_t type_bsize, index_t bound, index_t dim)
	: _node(node), _buffer(buffer), _type_bsize(type_bsize), _boundary_size(bound), _dim_size(dim)
	{
	}

	virtual void do_sync();

private:
	Node1D *_node;
	void *_buffer;
	index_t _type_bsize;
	index_t _boundary_size;
	index_t _dim_size;

	RemoteBoundary1D();
};

}
