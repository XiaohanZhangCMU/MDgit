#include <iostream>
#include <mpi.h>
#include "MachineInfo.h"

namespace StencilToolkit
{

MachineInfo &MachineInfo::get_instance()
{
	static MachineInfo info;
	return info;
}

void MachineInfo::error(const char *msg)
{
	std::cout << msg << std::endl;
	MPI_Abort(MPI_COMM_WORLD, -1);
}

MachineInfo::MachineInfo()
{
	int size, rank;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	_world_size = static_cast<index_t>(size);
	_world_rank = static_cast<index_t>(rank);
}

void MachineInfo::reduce(ReductionOp op, ReductionType type, void *buffer)
{
	MPI_Op mpi_op;
	switch (op)
	{
	case STK_SUM: mpi_op = MPI_SUM; break;
	case STK_MAX: mpi_op = MPI_MAX; break;
	}

	MPI_Datatype mpi_type;
	switch (type)
	{
	case STK_INT:    mpi_type = MPI_INT;    break;
	case STK_FLOAT:  mpi_type = MPI_FLOAT;  break;
	case STK_DOUBLE: mpi_type = MPI_DOUBLE; break;
	}

	MPI_Allreduce(MPI_IN_PLACE, buffer, 1, mpi_type, mpi_op, MPI_COMM_WORLD);
}

}
