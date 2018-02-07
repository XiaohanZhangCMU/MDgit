#include <cstring>
#include "LocalBoundary1D.h"

namespace StencilToolkit
{

void LocalBoundary1D::do_sync()
{
	char *c_buf = static_cast<char *>(_buffer);

	index_t dim_bsize = _dim_size * _type_bsize;
	index_t bound_bsize = _boundary_size * _type_bsize;

	std::memcpy(c_buf,
	            c_buf + dim_bsize,
	            bound_bsize);

	std::memcpy(c_buf + bound_bsize + dim_bsize,
	            c_buf + bound_bsize,
	            bound_bsize);
}

}
