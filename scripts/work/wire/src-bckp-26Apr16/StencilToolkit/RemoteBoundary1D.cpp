#include <mpi.h>
#include "RemoteBoundary1D.h"
#include <cstring>

namespace StencilToolkit
{

void RemoteBoundary1D::do_sync()
{
	if (_boundary_size == 0)
	{
		return;
	}

	char *c_buf = static_cast<char *>(_buffer);

	index_t dim_bsize = _dim_size * _type_bsize;
	index_t bound_bsize = _boundary_size * _type_bsize;

	// FIXME consider better implementation
	if (_node->size() == 1)
	{
		std::memcpy(c_buf,
		            c_buf + dim_bsize,
		            bound_bsize);

		std::memcpy(c_buf + bound_bsize + dim_bsize,
		            c_buf + bound_bsize,
		            bound_bsize);
	}
	else
	{
		int np1 = (*_node)(+1);
		int nm1 = (*_node)(-1);

		MPI_Request req[4];
		MPI_Status stat[4];

		MPI_Irecv(c_buf,                           bound_bsize, MPI_BYTE, nm1, 0, MPI_COMM_WORLD, &(req[0]));
		MPI_Isend(c_buf + dim_bsize,               bound_bsize, MPI_BYTE, np1, 0, MPI_COMM_WORLD, &(req[1]));
		MPI_Isend(c_buf + bound_bsize,             bound_bsize, MPI_BYTE, nm1, 1, MPI_COMM_WORLD, &(req[2]));
		MPI_Irecv(c_buf + bound_bsize + dim_bsize, bound_bsize, MPI_BYTE, np1, 1, MPI_COMM_WORLD, &(req[3]));
		MPI_Waitall(4, req, stat);
	}
}

}
