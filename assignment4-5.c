/***************************************************************************/
/* Template for Asssignment 4/5 ********************************************/
/* Team Names Here              **(*****************************************/
/***************************************************************************/

/* 

Compile with:
mpicc -g -Wall -o hw4 assignment4-5.c clcg4.c -lpthread
Run with:
mpirun -np 4 ./hw4

- Notes -

ticks are a user input. one tick is running through the whole grid once and 
applying the rules (including randomization if threshold says so)

Ghost rows are the rows representative of the boundaries of other ranks so 
your current rank can read them when determining the logic. 

pthreads seem to go right past mpi collectives like MPI_Barrier. using 
pthread barriers is a good way to do synchronisation. 

just call InitDefault(); in all threads, then call GenVal() with the row 
number to get the random value 0-1 for each cell for each row within each thread

I spoke with the professor in class today and he said that cell updates are 
based on neighbor cells that have been updated within the same tick.
Cell 0 updates -> Cell 1 updates based on updated_Cell_0
we do row by row, using previous changed values as new neighbor values.
    so cells updates are serial but row updates can be from previous tick

Randomness: random value for every cell each tick, if below some threshold, 
    apply a random state

*/

/***************************************************************************/
/* Includes ****************************************************************/
/***************************************************************************/

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<errno.h>
#include<math.h>

#include"clcg4.h"

#include<mpi.h>
#include<pthread.h>

// #define BGQ 1 // when running BG/Q, comment out when testing on mastiff

#ifdef BGQ
#include<hwi/include/bqc/A2_inlines.h>
#else
#define GetTimeBase MPI_Wtime            
#endif

/***************************************************************************/
/* Defines *****************************************************************/
/***************************************************************************/

#define ALIVE 1
#define DEAD  0

#define board_size 8//32768

/***************************************************************************/
/* Global Vars *************************************************************/
/***************************************************************************/

double g_time_in_secs = 0;
double g_processor_frequency = 1600000000.0; // processing speed for BG/Q
unsigned long long g_start_cycles=0;
unsigned long long g_end_cycles=0;

// You define these


/***************************************************************************/
/* Function Decs ***********************************************************/
/***************************************************************************/

// You define these
void allocate_mem(int*** arr, int h, int w);
void deallocate_mem(int*** arr, int h);
void print_board(int*** arr, int h, int w);


/***************************************************************************/
/* Function: Main **********************************************************/
/***************************************************************************/

int main(int argc, char *argv[])
{
//    int i = 0;
    int mpi_myrank;
    int mpi_commsize;
// Example MPI startup and using CLCG4 RNG
    MPI_Init( &argc, &argv);
    MPI_Comm_size( MPI_COMM_WORLD, &mpi_commsize);
    MPI_Comm_rank( MPI_COMM_WORLD, &mpi_myrank);
    
// Init 32,768 RNG streams - each rank has an independent stream
    InitDefault();
    
// Note, used the mpi_myrank to select which RNG stream to use.
// You must replace mpi_myrank with the right row being used.
// This just show you how to call the RNG.    
    // printf("Rank %d of %d has been started and a first Random Value of %lf\n", 
	   // mpi_myrank, mpi_commsize, GenVal(mpi_myrank));
    
    MPI_Barrier( MPI_COMM_WORLD );
    
// Insert your code
    // Following algorithm description from class
    // =========================================================
    // =========================================================
    // =========================================================
    // =========================================================
    // =========================================================
    int num_threads = 4;

    // if rank 0 / pthread0, start time with GetTimeBase() 
    if(mpi_myrank == 0){
        g_start_cycles = GetTimeBase();
        
        // sanity check
        printf("Board size: %d x %d\n", board_size, board_size);
        int cells_per_rank = board_size*(board_size/mpi_commsize);
        int rows_per_rank = board_size/mpi_commsize;
        printf("each of %d ranks will have %d cells, %d rows, %d ghost rows\n",
            mpi_commsize,cells_per_rank, rows_per_rank, 2);
        int cells_per_thread = cells_per_rank/num_threads;
        int rows_per_thread = rows_per_rank/num_threads;
        printf("each of %d threads will have %d cells, %d rows\n",
            num_threads,cells_per_thread, rows_per_thread);
    }

    // attempt to allocate dynamically in other functions (didnt work)
    // int h=2, w=2;
    // int ** board;
    // allocate_mem(&board, h, w);
    // print_board(&board, h, w);

    /*  1
        allocate my rank's chunk of the universe + space for ghost rows
            ghost rows are the rows at edges of rank boundaries
            allocating 2d array of the board for each rank
    */
    // allocate dynamically within main bc in functions was taking up time
    // declare main variables
    int num_ticks = 10;
    int rows_per_rank = board_size/mpi_commsize; 
    int rows_per_thread = rows_per_rank/num_threads;
    int **board;
    int *ghost_above;
    int *ghost_below;

    // allocate memory for the board
    board = (int **)malloc(rows_per_rank * sizeof(int *)); 
    for (int i=0; i<rows_per_rank; i++){
        board[i] = (int *)malloc(board_size * sizeof(int)); 
    }
    // allocate memory for the two ghost rows
    ghost_above = (int *)malloc(board_size * sizeof(int)); 
    ghost_below = (int *)malloc(board_size * sizeof(int)); 
    
    // initialize the board
    for(int i=0; i<rows_per_rank; i++){
        for(int j=0; j<board_size; j++){
            board[i][j] = ALIVE;
        }
    }
    // initialize the ghost row above
    for(int i=0; i<board_size; i++){
        ghost_above[i] = ALIVE;
        ghost_below[i] = ALIVE;
    }
  

    // // print the board for sanity
    // for(int i=0; i<rows_per_rank; i++){
    //     for(int j=0; j<board_size; j++){ 
    //         printf("%d ", board[i][j]); 
    //     }
    //     printf("\n");
    // }
    // printf("\n");
    // // print the ghost row above
    // for(int i=0; i<board_size; i++){ 
    //     printf("%d ", ghost_above[i]); 
    // }
    // printf("\n");
    // // print the ghost row below
    // for(int i=0; i<board_size; i++){ 
    //     printf("%d ", ghost_below[i]); 
    // }
    // printf("\n");

    /*  2
        Create Pthreads here. All threads should go into for loop
    */





    /*  3
        For all number of ticks, complete a round of the GOL
    */
    for(int tick=0; tick<num_ticks; tick++){

        /*  4
            Exchange row data with MPI ranks 
            using MPI_Isend/Irecv from thread 0 w/i each MPI rank.
            Yes, you must correctly MPI_Test or Wait to make sure
            messages operations correctly complete.

            [Note: have only 1 MPI rank/pthread perform ALL MPI
            operations per rank/thread group. Dont’ allow multiple
            threads to perform any MPI operations within MPI
            rank/thread group.]

        */  



        /*  5
            HERE each PTHREAD can process a row:
            - update universe making sure to use the
                correct row RNG stream
            - factor in Threshold percentage as described
            - use the right "ghost" row data at rank boundaries
            - keep track of total number of ALIVE cells per tick
                across all threads w/i a MPI rank group.
            - use pthread_mutex_trylock around shared counter
                variables **if needed**.
        */
    }

    /*  6
        MPI_Reduce( Sum all ALIVE Cells Counts For Each Tick);
        - Here, you will have vector of 256 ALIVE cell sum
        values which is the total number of ALIVE cells
        at each tick, t for all 256 ticks.
    */


    /*  7
        if rank 0 / thread 0, end time with GetTimeBase().
        if needed by experiment,
        perform output of 32Kx32K cell universe using MPI_file_write_at;
        collect I/O performance stats using GetTimeBase() from
        rank 0 / thread 0;
        if needed by experiment,
        construct 1Kx1K heatmap of 32Kx32K cell universe using MPI
        collective operations of your choice. Rank 0 will output
        the heatmap to a standard Unix file given it’s small 1 to 4MB size.
        Make sure you an import data for graphing.
        if rank 0, print ALIVE tick stats and compute (I/O if needed)
        performance stats.
    */







    // Deallocate memory
    
    for (int i = 0; i<rows_per_rank; i++){
        free(board[i]);
    }
    free(board); 
    free(ghost_above);
    free(ghost_below);

    
    

 
    // =========================================================
    // =========================================================
    // =========================================================
    // =========================================================
    // =========================================================
// END -Perform a barrier and then leave MPI
    MPI_Barrier( MPI_COMM_WORLD );
    MPI_Finalize();
    return 0;
}

/***************************************************************************/
/* Other Functions - You write as part of the assignment********************/
/***************************************************************************/

void allocate_mem(int*** arr, int h, int w){
    *arr = (int**)malloc(h*sizeof(int*));
    for(int i=0; i<h; i++){
        (*arr)[i] = (int*)malloc(w*sizeof(int));
    }
} 

void deallocate_mem(int*** arr, int h){
    for (int i = 0; i<h; i++){
        free((*arr)[i]);
    }
    free(*arr); 
}

void initialize_array(int** arr, int h, int w){
    for(int i=0; i<h; i++){
        for(int j=0; j<w; j++){
            arr[i][j] = 0;
        }
    }
}

void print_board(int*** arr, int h, int w){
    for(int i=0; i<h; i++){
        for(int j=0; j<w; j++){
            printf("%d", (*arr)[i][j]);
        }
        printf("\n");
    }
}