#pragma once

#include <cstring>
#include "Index2D.hpp"
#include "Uncopyable.hpp"

namespace StencilToolkit
{

template <typename T>
class Grid2D : private Uncopyable
{
public:
	Grid2D(Index2D &idx,
	       index_t bound_x = 0,
	       index_t bound_y = 0)
	: _idx(idx)
	{
		init(bound_x, bound_y);
	}

	~Grid2D()
	{
		delete[] _buffer;
		_buffer = nullptr;

		delete _boundary;
		_boundary = nullptr;
	}

	T &operator()(index_t x, index_t y)
	{
		return _buffer[calc_index(x, y)];
	}

	void set_boundary_periodic()
	{
		_boundary->do_sync();
	}

private:
	void init(index_t bound_x, index_t bound_y)
	{
		index_t size = ((_idx.get_range_x().end() - _idx.get_range_x().begin()) + (2 * bound_x))
		             * ((_idx.get_range_y().end() - _idx.get_range_y().begin()) + (2 * bound_y));
		_buffer = new T[size];
		std::memset(_buffer, 0, size * sizeof(T));

		_boundary = _idx.alloc_boundary(_buffer, sizeof(T), bound_x, bound_y);

		_diff_y = (_idx.get_range_x().end() - _idx.get_range_x().begin()) + (2 * bound_x);

		_G2L_diff = (bound_x - _idx.get_range_x().begin()) +
		           ((bound_y - _idx.get_range_y().begin()) * _diff_y);
	}

	index_t calc_index(index_t x, index_t y)
	{
		return _G2L_diff + x + (y * _diff_y);
	}

	T *_buffer;

	Index2D &_idx;

	Boundary *_boundary;

	index_t _diff_y;

	index_t _G2L_diff;

	Grid2D();
};

}
