#include "partition/MPCNS_Pre_Datatrans.h"

//=================================================================================================
// MPCNS_IF_MPI=1才会生效
#if MPCNS_IF_MPI == 1

namespace DATATRANS
{
    void mpi_initial(int arg, char **argv)
    {
        MPI_Init(&arg, &argv);
    }
    void mpi_finalize()
    {
        MPI_Finalize();
    }
    void mpi_barrier()
    {
        MPI_Barrier(MPI_COMM_WORLD);
    }
    void mpi_data_send(int tar_id, int send_flag, double *data, int length, MPI_Request *status)
    {
        MPI_Isend(data, length, MPI_DOUBLE, tar_id, send_flag, MPI_COMM_WORLD, status);
    }
    void mpi_data_recv(int tar_id, int recv_flag, double *data, int length, MPI_Request *status)
    {
        MPI_Irecv(data, length, MPI_DOUBLE, tar_id, recv_flag, MPI_COMM_WORLD, status);
    }
    void mpi_rank(int *rank)
    {
        MPI_Comm_rank(MPI_COMM_WORLD, rank);
    }
    void mpi_size(int *size)
    {
        MPI_Comm_size(MPI_COMM_WORLD, size);
    }
    void mpi_wait(int &count, MPI_Request *request_s, MPI_Status *status_s)
    {
        MPI_Waitall(count, request_s, status_s);
    }
}
#endif
//=================================================================================================
