#include <mpi.h>
#include "RemoteBoundary3D.h"
#include <cstring>

namespace StencilToolkit
{

void RemoteBoundary3D::do_sync()
{
	char *c_buf = static_cast<char *>(_buffer);
	MPI_Datatype *types = static_cast<MPI_Datatype *>(_boundary_type);

	index_t diff_y = _dim_size_x + (2 * _boundary_size_x);
	index_t diff_z = diff_y * (_dim_size_y + (2 * _boundary_size_y));
	index_t start_index = _boundary_size_x + (_boundary_size_y * diff_y) + (_boundary_size_z * diff_z);

	// FIXME more efficient implementation
	if (_boundary_size_x != 0)
	{
		if (_node->size_x() == 1)
		{
			for (int i = 0; i < _dim_size_z; i++)
			{
				for (int j = 0; j < _dim_size_y; j++)
				{
					std::memcpy(c_buf + calc_bpos(_type_bsize, -_boundary_size_x,              j, i, start_index, diff_y, diff_z),
					            c_buf + calc_bpos(_type_bsize, _dim_size_x - _boundary_size_x, j, i, start_index, diff_y, diff_z),
					            _boundary_size_x * _type_bsize);

					std::memcpy(c_buf + calc_bpos(_type_bsize, _dim_size_x, j, i, start_index, diff_y, diff_z),
					            c_buf + calc_bpos(_type_bsize, 0,           j, i, start_index, diff_y, diff_z),
					            _boundary_size_x * _type_bsize);
				}
			}
		}
		else
		{
			int np1_x = (*_node)(+1,  0,  0);
			int nm1_x = (*_node)(-1,  0,  0);

			for (int i = 0; i < _dim_size_z; i++)
			{
				MPI_Request req[4];
				MPI_Irecv(c_buf + calc_bpos(_type_bsize, -_boundary_size_x,              0, i, start_index, diff_y, diff_z), 1, types[0], nm1_x, 0, MPI_COMM_WORLD, &(req[0]));
				MPI_Isend(c_buf + calc_bpos(_type_bsize, _dim_size_x - _boundary_size_x, 0, i, start_index, diff_y, diff_z), 1, types[0], np1_x, 0, MPI_COMM_WORLD, &(req[1]));
				MPI_Isend(c_buf + calc_bpos(_type_bsize, 0,                              0, i, start_index, diff_y, diff_z), 1, types[0], nm1_x, 1, MPI_COMM_WORLD, &(req[2]));
				MPI_Irecv(c_buf + calc_bpos(_type_bsize, _dim_size_x,                    0, i, start_index, diff_y, diff_z), 1, types[0], np1_x, 1, MPI_COMM_WORLD, &(req[3]));

				MPI_Status stat[4];
				MPI_Waitall(4, req, stat);
			}
		}
	}

	if (_boundary_size_y != 0)
	{
		if (_node->size_y() == 1)
		{
			index_t alloc_size_x = _dim_size_x + (2 * _boundary_size_x);
			index_t alloc_size_y = _dim_size_y + (2 * _boundary_size_y);

			for (int i = 0; i < _dim_size_z; i++)
			{
				std::memcpy(c_buf + calc_bpos(_type_bsize, -_boundary_size_x, -_boundary_size_y,             i, start_index, diff_y, diff_z),
				            c_buf + calc_bpos(_type_bsize, -_boundary_size_x, _dim_size_y -_boundary_size_y, i, start_index, diff_y, diff_z),
				            alloc_size_x * _boundary_size_y * _type_bsize);

				std::memcpy(c_buf + calc_bpos(_type_bsize, -_boundary_size_x, _dim_size_y, i, start_index, diff_y, diff_z),
				            c_buf + calc_bpos(_type_bsize, -_boundary_size_x, 0,           i, start_index, diff_y, diff_z),
				            alloc_size_x * _boundary_size_y * _type_bsize);
			}
		}
		else
		{
			int np1_y = (*_node)( 0, +1,  0);
			int nm1_y = (*_node)( 0, -1,  0);

			MPI_Request req[4];
			MPI_Irecv(c_buf + calc_bpos(_type_bsize, -_boundary_size_x, -_boundary_size_y,             0, start_index, diff_y, diff_z), 1, types[1], nm1_y, 0, MPI_COMM_WORLD, &(req[0]));
			MPI_Isend(c_buf + calc_bpos(_type_bsize, -_boundary_size_x, _dim_size_y -_boundary_size_y, 0, start_index, diff_y, diff_z), 1, types[1], np1_y, 0, MPI_COMM_WORLD, &(req[1]));
			MPI_Isend(c_buf + calc_bpos(_type_bsize, -_boundary_size_x, 0,                             0, start_index, diff_y, diff_z), 1, types[1], nm1_y, 1, MPI_COMM_WORLD, &(req[2]));
			MPI_Irecv(c_buf + calc_bpos(_type_bsize, -_boundary_size_x, _dim_size_y,                   0, start_index, diff_y, diff_z), 1, types[1], np1_y, 1, MPI_COMM_WORLD, &(req[3]));

			MPI_Status stat[4];
			MPI_Waitall(4, req, stat);
		}
	}

	if (_boundary_size_z != 0)
	{
		if (_node->size_z() == 1)
		{
			index_t alloc_size_x = _dim_size_x + (2 * _boundary_size_x);
			index_t alloc_size_y = _dim_size_y + (2 * _boundary_size_y);

			std::memcpy(c_buf + calc_bpos(_type_bsize, -_boundary_size_x, -_boundary_size_y, -_boundary_size_z,             start_index, diff_y, diff_z),
			            c_buf + calc_bpos(_type_bsize, -_boundary_size_x, -_boundary_size_y, _dim_size_z -_boundary_size_z, start_index, diff_y, diff_z),
			            alloc_size_x * alloc_size_y * _boundary_size_z * _type_bsize);

			std::memcpy(c_buf + calc_bpos(_type_bsize, -_boundary_size_x, -_boundary_size_y, _dim_size_z, start_index, diff_y, diff_z),
			            c_buf + calc_bpos(_type_bsize, -_boundary_size_x, -_boundary_size_y, 0,           start_index, diff_y, diff_z),
			            alloc_size_x * alloc_size_y * _boundary_size_z * _type_bsize);
		}
		else
		{
			int np1_z = (*_node)( 0,  0, +1);
			int nm1_z = (*_node)( 0,  0, -1);

			MPI_Request req[4];
			MPI_Irecv(c_buf + calc_bpos(_type_bsize, -_boundary_size_x, -_boundary_size_y, -_boundary_size_z,             start_index, diff_y, diff_z),
			                            1, types[2], nm1_z, 0, MPI_COMM_WORLD, &(req[0]));
			MPI_Isend(c_buf + calc_bpos(_type_bsize, -_boundary_size_x, -_boundary_size_y, _dim_size_z -_boundary_size_z, start_index, diff_y, diff_z),
			                            1, types[2], np1_z, 0, MPI_COMM_WORLD, &(req[1]));
			MPI_Isend(c_buf + calc_bpos(_type_bsize, -_boundary_size_x, -_boundary_size_y, 0,                             start_index, diff_y, diff_z),
			                            1, types[2], nm1_z, 1, MPI_COMM_WORLD, &(req[2]));
			MPI_Irecv(c_buf + calc_bpos(_type_bsize, -_boundary_size_x, -_boundary_size_y, _dim_size_z,                   start_index, diff_y, diff_z),
			                            1, types[2], np1_z, 1, MPI_COMM_WORLD, &(req[3]));

			MPI_Status stat[4];
			MPI_Waitall(4, req, stat);
		}
	}
}

void RemoteBoundary3D::init_type()
{
	MPI_Datatype *types = new MPI_Datatype[3];

	index_t diff_y = _dim_size_x + (2 * _boundary_size_x);
	index_t diff_z = diff_y * (_dim_size_y + (2 * _boundary_size_y));
	index_t alloc_size_x = _dim_size_x + (2 * _boundary_size_x);
	index_t alloc_size_y = _dim_size_y + (2 * _boundary_size_y);

	MPI_Type_vector(_dim_size_y, _boundary_size_x * _type_bsize, diff_y * _type_bsize, MPI_BYTE, &(types[0]));
	MPI_Type_commit(&(types[0]));
// FIXME declare 3D vector DDT
//	MPI_Type_vector(_dim_size_z, 1, diff_z * _type_bsize, temp_type, &(types[0]));
//	MPI_Type_commit(&(types[0]));

	MPI_Type_vector(_dim_size_z, alloc_size_x * _boundary_size_y * _type_bsize, diff_z * _type_bsize, MPI_BYTE, &(types[1]));
	MPI_Type_commit(&(types[1]));

	MPI_Type_contiguous(alloc_size_x * alloc_size_y * _boundary_size_z * _type_bsize, MPI_BYTE, &(types[2]));
	MPI_Type_commit(&(types[2]));

	_boundary_type = types;
}

void RemoteBoundary3D::finalize_type()
{
	MPI_Datatype *types = static_cast<MPI_Datatype *>(_boundary_type);

	MPI_Type_free(&(types[0]));
	MPI_Type_free(&(types[1]));
	MPI_Type_free(&(types[2]));

	delete[] types;
}

}
