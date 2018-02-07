#pragma once

#include "Range.hpp"

namespace StencilToolkit
{

class SerialRange : public Range
{
public:
	SerialRange(index_t size)
	: _size(size)
	{
	}

	virtual index_t begin()
	{
		return 0;
	}

	virtual index_t end()
	{
		return _size;
	}

private:
	index_t _size;

	SerialRange();
};

}
