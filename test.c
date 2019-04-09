/* This is an interactive version of cpi */
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


int main(int argc,char *argv[])
{

    int  num_ranks, my_rank;
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
    MPI_Comm_size(MPI_COMM_WORLD,&num_ranks);    
    MPI_Status status;

    MPI_File fh;


    if (my_rank == 0) {
        MPI_File_open(MPI_COMM_SELF, "test.txt",MPI_MODE_CREATE | MPI_MODE_WRONLY,MPI_INFO_NULL,&fh);
        
        int length = 10;
        long* y = malloc(length*sizeof(long));
        for(int i=0;i<length;i++){
            y[i] = i*i;
            MPI_File_write(fh, y[i], 1, MPI_LONG, &status);
        }
        
        
        //        fclose(f);
        MPI_File_close(&fh);
    }
    else {
        // do nothing
    }


    MPI_Finalize();
    return 0;
}