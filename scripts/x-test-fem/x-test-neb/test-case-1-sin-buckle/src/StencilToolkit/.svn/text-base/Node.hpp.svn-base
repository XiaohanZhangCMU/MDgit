#pragma once

#include "common.h"
#include "MachineInfo.h"
#include "reduce.h"
#include "Uncopyable.hpp"

namespace StencilToolkit
{

class Node : private Uncopyable
{
public:
	virtual bool is_master()
	{
		return true;
	}

	// FIXME consider better interface, implementation
	void reduce(ReductionOp op, ReductionType type, void *buffer)
	{
		MachineInfo::get_instance().reduce(op, type, buffer);
	}
};

}
