#pragma once

#include "common.h"
#include "Boundary.hpp"

namespace StencilToolkit
{

class LocalBoundary1D : public Boundary
{
public:
	LocalBoundary1D(void *buffer, index_t type_bsize, index_t bound, index_t dim)
	: _buffer(buffer), _type_bsize(type_bsize), _boundary_size(bound), _dim_size(dim)
	{
	}

	virtual void do_sync();

private:
	void *_buffer;
	index_t _type_bsize;
	index_t _boundary_size;
	index_t _dim_size;

	LocalBoundary1D();
};

}
