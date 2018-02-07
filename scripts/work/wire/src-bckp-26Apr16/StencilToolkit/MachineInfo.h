#pragma once

#include "common.h"
#include "reduce.h"

namespace StencilToolkit
{

class MachineInfo
{
public:
	static MachineInfo &get_instance();
	static void error(const char *msg);

	index_t get_size()
	{
		return _world_size;
	}

	index_t get_rank()
	{
		return _world_rank;
	}

	void reduce(ReductionOp op, ReductionType type, void *buffer);

private:
	MachineInfo();
	MachineInfo(const MachineInfo &);
	MachineInfo &operator=(const MachineInfo &);

	index_t _world_size;
	index_t _world_rank;
};

}
