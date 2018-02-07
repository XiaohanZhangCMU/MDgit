#pragma once

#include "MachineInfo.h"
#include "Range.hpp"

namespace StencilToolkit
{

class BlockRange : public Range
{
public:
	BlockRange(index_t dim_size, index_t node_size, index_t node_rank)
	{
		index_t chunk_size = dim_size / node_size;
		if (chunk_size == 0)
		{
			MachineInfo::error("[error] too small size");
		}

		_begin = chunk_size * node_rank;
		if (node_rank == (node_size - 1))
		{
			_end = dim_size;
		}
		else
		{
			_end = _begin + chunk_size;
		}
	}

	virtual index_t begin()
	{
		return _begin;
	}

	virtual index_t end()
	{
		return _end;
	}

private:
	index_t _begin;
	index_t _end;

	BlockRange();
};

}
