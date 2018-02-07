#include <cstring>
#include "LocalBoundary2D.h"

namespace StencilToolkit
{

void LocalBoundary2D::do_sync()
{
	char *c_buf = static_cast<char *>(_buffer);

	index_t diff_y = _dim_size_x + (2 * _boundary_size_x);
	index_t start_index = _boundary_size_x + (_boundary_size_y * diff_y);

	for (int i = 0; i < _dim_size_y; i++)
	{
		std::memcpy(c_buf + calc_bpos(_type_bsize, -_boundary_size_x,              i, start_index, diff_y),
		            c_buf + calc_bpos(_type_bsize, _dim_size_x - _boundary_size_x, i, start_index, diff_y),
		            _boundary_size_x * _type_bsize);

		std::memcpy(c_buf + calc_bpos(_type_bsize, _dim_size_x, i, start_index, diff_y),
		            c_buf + calc_bpos(_type_bsize, 0,           i, start_index, diff_y),
		            _boundary_size_x * _type_bsize);
	}

	std::memcpy(c_buf + calc_bpos(_type_bsize, -_boundary_size_x, -_boundary_size_y,             start_index, diff_y),
	            c_buf + calc_bpos(_type_bsize, -_boundary_size_x, _dim_size_y -_boundary_size_y, start_index, diff_y),
	            (_dim_size_x + (2 * _boundary_size_x)) * _boundary_size_y * _type_bsize);

	std::memcpy(c_buf + calc_bpos(_type_bsize, -_boundary_size_x, _dim_size_y, start_index, diff_y),
	            c_buf + calc_bpos(_type_bsize, -_boundary_size_x, 0,           start_index, diff_y),
	            (_dim_size_x + (2 * _boundary_size_x)) * _boundary_size_y * _type_bsize);
}

}
