#ifndef ICS_MPI_H
#define ICS_MPI_H
#include <vector>

#include "mpi.h"
#include "vector3.hpp"

void mpi_sum_Vector3(Vector3* vs_in, Vector3* vs_out, numeric::size_type* len, MPI_Datatype* datatype)
{
    for(numeric::size_type i=0; i < *len; ++i)
    {
        vs_out[i] += vs_in[i]; 
    }
}


class Communicator
{
public:
    typedef int mpi_int_type; 

public:
    Communicator()
    {
        MPI_Init(NULL, NULL);
        world = MPI_COMM_WORLD;
        MPI_Comm_size(world, &world_size);
        MPI_Comm_rank(world, &world_rank);

        time_start = MPI_Wtime();
        
        regist_Vector3();
    }

    ~Communicator()
    {
        degregist_Vector3();

        MPI_Finalize();
    }

    mpi_int_type rank()
    {
        return world_rank;
    }
    mpi_int_type size()
    {
        return world_size;
    }

    MPI_Comm& getworld()
    {
        return world;
    }
    void barrier()
    {
        MPI_Barrier(world);
    }

    void broadcast(numeric::value_type& para)
    {
        MPI_Bcast(&para, 1, MPI_DOUBLE, 0, world);
    }
    
    void find_max(const numeric::value_type& data, numeric::value_type& data_max)
    {
        MPI_Allreduce(&data, &data_max, 1, MPI_DOUBLE, MPI_MAX,
              world);
    }

    void find_min(const numeric::value_type& data, numeric::value_type& data_min)
    {
        MPI_Allreduce(&data, &data_min, 1, MPI_DOUBLE, MPI_MIN,
              world);
    }

    void scatter_Vector3s(const std::vector<Vector3>& data_global, std::vector<Vector3>& data_local, std::vector<numeric::size_type>& counts, std::vector<numeric::size_type>& displs, mpi_int_type from_rank=0)
    {
        MPI_Scatterv(data_global.data(), counts.data(), displs.data(), mpi_Vector3, data_local.data(), data_local.size(), mpi_Vector3, from_rank, world);
    }
    
    void allgather(const mpi_int_type& data, std::vector<mpi_int_type>& data_list)
    {
        MPI_Allgather(&data, 1, MPI_INT, data_list.data(), 1, MPI_INT, world);
    }

    void sum_Vector3s(std::vector<Vector3>& input, std::vector<Vector3>& output, mpi_int_type to_rank=0)
    {
        MPI_Reduce(input.data(), output.data(), input.size(), mpi_Vector3, mpi_Vector3_sum, to_rank, world);
    }
    
    void dump_excution_time()
    {
        if(world_rank == 0)
        {
            std::cout << "Excution Time: ";
            std::cout << (MPI_Wtime() - time_start) << "\n";
        }
    }

private:
    void regist_Vector3()
    {
        MPI_Type_contiguous(3,MPI_DOUBLE,&mpi_Vector3);
        MPI_Type_commit(&mpi_Vector3);

        MPI_Op_create(mpi_sum_Vector3, true, &mpi_Vector3_sum);
    }
    void degregist_Vector3()
    {
        MPI_Op_free(&mpi_Vector3_sum);
        MPI_Type_free(&mpi_Vector3);
    }

private:
    MPI_Comm world;
    mpi_int_type world_size;
    mpi_int_type world_rank;

    MPI_Datatype mpi_Vector3;
    MPI_Op mpi_Vector3_sum;

    double time_start;
};


#endif
