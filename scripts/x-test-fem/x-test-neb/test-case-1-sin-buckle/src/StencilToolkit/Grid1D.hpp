#pragma once

#include <cstring>
#include "Boundary.hpp"
#include "Index1D.hpp"
#include "Uncopyable.hpp"

namespace StencilToolkit
{

template <typename T>
class Grid1D : private Uncopyable
{
public:
	Grid1D(Index1D &idx, index_t bound = 0)
	: _idx(idx)
	{
		init(bound);
	}

	~Grid1D()
	{
		delete[] _buffer;
		_buffer = nullptr;

		delete _boundary;
		_boundary = nullptr;
	}

	T &operator()(index_t x)
	{
		return _buffer[_G2L_diff + x];
	}

	void set_boundary_periodic()
	{
		_boundary->do_sync();
	}

private:
	void init(index_t bound)
	{
		index_t size = ((_idx.get_range().end() - _idx.get_range().begin()) + (2 * bound));
		_buffer = new T[size];
		std::memset(_buffer, 0, size * sizeof(T));

		_boundary = _idx.alloc_boundary(_buffer, sizeof(T), bound);

		_G2L_diff = bound - _idx.get_range().begin();
	}

	T *_buffer;

	Index1D &_idx;

	Boundary *_boundary;

	index_t _G2L_diff;

	Grid1D();
};

}
