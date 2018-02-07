#include <cstring>
#include "LocalBoundary3D.h"

namespace StencilToolkit
{

void LocalBoundary3D::do_sync()
{
	char *c_buf = static_cast<char *>(_buffer);

	index_t diff_y = _dim_size_x + (2 * _boundary_size_x);
	index_t diff_z = diff_y * (_dim_size_y + (2 * _boundary_size_y));
	index_t start_index = _boundary_size_x + (_boundary_size_y * diff_y) + (_boundary_size_z * diff_z);

	for (int i = 0; i < _dim_size_z; i++)
	{
		for (int j = 0; j < _dim_size_y; j++)
		{
			std::memcpy(c_buf + calc_bpos(_type_bsize, -_boundary_size_x,              j, i, start_index, diff_y, diff_z),
			            c_buf + calc_bpos(_type_bsize, _dim_size_x - _boundary_size_x, j, i, start_index, diff_y, diff_z),
			            _boundary_size_x * _type_bsize);

			std::memcpy(c_buf + calc_bpos(_type_bsize, _dim_size_x, j, i, start_index, diff_y, diff_z),
			            c_buf + calc_bpos(_type_bsize, 0,           j, i, start_index, diff_y, diff_z),
			            _boundary_size_x * _type_bsize);
		}
	}

	index_t alloc_size_x = _dim_size_x + (2 * _boundary_size_x);
	index_t alloc_size_y = _dim_size_y + (2 * _boundary_size_y);

	for (int i = 0; i < _dim_size_z; i++)
	{
		std::memcpy(c_buf + calc_bpos(_type_bsize, -_boundary_size_x, -_boundary_size_y,             i, start_index, diff_y, diff_z),
		            c_buf + calc_bpos(_type_bsize, -_boundary_size_x, _dim_size_y -_boundary_size_y, i, start_index, diff_y, diff_z),
		            alloc_size_x * _boundary_size_y * _type_bsize);

		std::memcpy(c_buf + calc_bpos(_type_bsize, -_boundary_size_x, _dim_size_y, i, start_index, diff_y, diff_z),
		            c_buf + calc_bpos(_type_bsize, -_boundary_size_x, 0,           i, start_index, diff_y, diff_z),
		            alloc_size_x * _boundary_size_y * _type_bsize);
	}

	std::memcpy(c_buf + calc_bpos(_type_bsize, -_boundary_size_x, -_boundary_size_y, -_boundary_size_z,             start_index, diff_y, diff_z),
	            c_buf + calc_bpos(_type_bsize, -_boundary_size_x, -_boundary_size_y, _dim_size_z -_boundary_size_z, start_index, diff_y, diff_z),
	            alloc_size_x * alloc_size_y * _boundary_size_z * _type_bsize);

	std::memcpy(c_buf + calc_bpos(_type_bsize, -_boundary_size_x, -_boundary_size_y, _dim_size_z, start_index, diff_y, diff_z),
	            c_buf + calc_bpos(_type_bsize, -_boundary_size_x, -_boundary_size_y, 0,           start_index, diff_y, diff_z),
	            alloc_size_x * alloc_size_y * _boundary_size_z * _type_bsize);
}

}
