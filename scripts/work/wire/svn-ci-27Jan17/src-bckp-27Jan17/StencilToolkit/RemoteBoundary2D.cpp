#include <mpi.h>
#include "RemoteBoundary2D.h"
#include <cstring>

namespace StencilToolkit
{

void RemoteBoundary2D::do_sync()
{
	char *c_buf = static_cast<char *>(_buffer);
	MPI_Datatype *types = static_cast<MPI_Datatype *>(_boundary_type);

	index_t diff_y = _dim_size_x + (2 * _boundary_size_x);
	index_t start_index = _boundary_size_x + (_boundary_size_y * diff_y);

	if (_boundary_size_x != 0)
	{
		if (_node->size_x() == 1)
		{
			for (int i = 0; i < _dim_size_y; i++)
			{
				std::memcpy(c_buf + calc_bpos(_type_bsize, -_boundary_size_x,              i, start_index, diff_y),
				            c_buf + calc_bpos(_type_bsize, _dim_size_x - _boundary_size_x, i, start_index, diff_y),
				            _boundary_size_x * _type_bsize);

				std::memcpy(c_buf + calc_bpos(_type_bsize, _dim_size_x, i, start_index, diff_y),
				            c_buf + calc_bpos(_type_bsize, 0,           i, start_index, diff_y),
				            _boundary_size_x * _type_bsize);
			}
		}
		else
		{
			int np1_x = (*_node)(+1,  0);
			int nm1_x = (*_node)(-1,  0);

			MPI_Request req[4];
			MPI_Irecv(c_buf + calc_bpos(_type_bsize, -_boundary_size_x,              0, start_index, diff_y), 1, types[0], nm1_x, 0, MPI_COMM_WORLD, &(req[0]));
			MPI_Isend(c_buf + calc_bpos(_type_bsize, _dim_size_x - _boundary_size_x, 0, start_index, diff_y), 1, types[0], np1_x, 0, MPI_COMM_WORLD, &(req[1]));
			MPI_Isend(c_buf + calc_bpos(_type_bsize, 0,                              0, start_index, diff_y), 1, types[0], nm1_x, 1, MPI_COMM_WORLD, &(req[2]));
			MPI_Irecv(c_buf + calc_bpos(_type_bsize, _dim_size_x,                    0, start_index, diff_y), 1, types[0], np1_x, 1, MPI_COMM_WORLD, &(req[3]));

			MPI_Status stat[4];
			MPI_Waitall(4, req, stat);
		}
	}

	if (_boundary_size_y != 0)
	{
		if (_node->size_y() == 1)
		{
			std::memcpy(c_buf + calc_bpos(_type_bsize, -_boundary_size_x, -_boundary_size_y,             start_index, diff_y),
			            c_buf + calc_bpos(_type_bsize, -_boundary_size_x, _dim_size_y -_boundary_size_y, start_index, diff_y),
			            (_dim_size_x + (2 * _boundary_size_x)) * _boundary_size_y * _type_bsize);

			std::memcpy(c_buf + calc_bpos(_type_bsize, -_boundary_size_x, _dim_size_y, start_index, diff_y),
			            c_buf + calc_bpos(_type_bsize, -_boundary_size_x, 0,           start_index, diff_y),
			            (_dim_size_x + (2 * _boundary_size_x)) * _boundary_size_y * _type_bsize);
		}
		else
		{
			int np1_y = (*_node)( 0, +1);
			int nm1_y = (*_node)( 0, -1);

			MPI_Request req[4];
			MPI_Irecv(c_buf + calc_bpos(_type_bsize, -_boundary_size_x, -_boundary_size_y,             start_index, diff_y), 1, types[1], nm1_y, 0, MPI_COMM_WORLD, &(req[0]));
			MPI_Isend(c_buf + calc_bpos(_type_bsize, -_boundary_size_x, _dim_size_y -_boundary_size_y, start_index, diff_y), 1, types[1], np1_y, 0, MPI_COMM_WORLD, &(req[1]));
			MPI_Isend(c_buf + calc_bpos(_type_bsize, -_boundary_size_x, 0,                             start_index, diff_y), 1, types[1], nm1_y, 1, MPI_COMM_WORLD, &(req[2]));
			MPI_Irecv(c_buf + calc_bpos(_type_bsize, -_boundary_size_x, _dim_size_y,                   start_index, diff_y), 1, types[1], np1_y, 1, MPI_COMM_WORLD, &(req[3]));

			MPI_Status stat[4];
			MPI_Waitall(4, req, stat);
		}
	}
}

void RemoteBoundary2D::init_type()
{
	index_t diff_y = _dim_size_x + (2 * _boundary_size_x);
	
	MPI_Datatype *types = new MPI_Datatype[2];

	MPI_Type_vector(_dim_size_y, _boundary_size_x * _type_bsize, diff_y * _type_bsize, MPI_BYTE, &(types[0]));
	MPI_Type_commit(&(types[0]));

	MPI_Type_contiguous(diff_y * _boundary_size_y * _type_bsize, MPI_BYTE, &(types[1]));
	MPI_Type_commit(&(types[1]));

	_boundary_type = types;
}

void RemoteBoundary2D::finalize_type()
{
	MPI_Datatype *types = static_cast<MPI_Datatype *>(_boundary_type);

	MPI_Type_free(&(types[0]));
	MPI_Type_free(&(types[1]));

	delete[] types;
}

}
