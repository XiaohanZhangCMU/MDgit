#pragma once

#include <cstring>
#include "Index3D.hpp"
#include "Uncopyable.hpp"

namespace StencilToolkit
{

template <typename T>
class Grid3D : private Uncopyable
{
public:
	Grid3D(Index3D &idx,
	       index_t bound_x = 0,
	       index_t bound_y = 0,
	       index_t bound_z = 0)
	: _idx(idx)
	{
		init(bound_x, bound_y, bound_z);
	}

	~Grid3D()
	{
		delete[] _buffer;
		_buffer = nullptr;

		delete _boundary;
		_boundary = nullptr;
	}

	T &operator()(index_t x, index_t y, index_t z)
	{
		return _buffer[calc_index(x, y, z)];
	}

	void set_boundary_periodic()
	{
		_boundary->do_sync();
	}

private:
	void init(index_t bound_x, index_t bound_y, index_t bound_z)
	{
		index_t size = ((_idx.get_range_x().end() - _idx.get_range_x().begin()) + (2 * bound_x))
		             * ((_idx.get_range_y().end() - _idx.get_range_y().begin()) + (2 * bound_y))
		             * ((_idx.get_range_z().end() - _idx.get_range_z().begin()) + (2 * bound_z));
		_buffer = new T[size];
		std::memset(_buffer, 0, size * sizeof(T));

		_boundary = _idx.alloc_boundary(_buffer, sizeof(T), bound_x, bound_y, bound_z);

		index_t dim_size_x = _idx.get_range_x().end() - _idx.get_range_x().begin();
		index_t dim_size_y = _idx.get_range_y().end() - _idx.get_range_y().begin();

		_diff_y = dim_size_x + (2 * bound_x);
		_diff_z = _diff_y * (dim_size_y + (2 * bound_y));

		_G2L_diff = (bound_x - _idx.get_range_x().begin()) +
		           ((bound_y - _idx.get_range_y().begin()) * _diff_y) +
		           ((bound_z - _idx.get_range_z().begin()) * _diff_z);
	}

	index_t calc_index(index_t x, index_t y, index_t z)
	{
		return _G2L_diff + x + (y * _diff_y) + (z * _diff_z);
	}

	T *_buffer;

	Index3D &_idx;

	Boundary *_boundary;

	index_t _diff_y;
	index_t _diff_z;

	index_t _G2L_diff;

	Grid3D();
};

}
