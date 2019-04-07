/***************************************************************************/
/* Template for Asssignment 4/5 ********************************************/
/* Austin Egri, Henry Grover    ********************************************/
/***************************************************************************/

/* 

Compile with:
mpicc -g -Wall -o hw4.out assignment4-5.c clcg4.c -lpthread
Run with:
mpirun -np 4 ./hw4.out

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

#define board_size 16//32768

/***************************************************************************/
/* Global Vars *************************************************************/
/***************************************************************************/

/*********************************/
// // MUTEX - ING
// pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;

// /** START MUTEX HERE  **/
// pthread_mutex_lock( &mutex );

// /* THIS CODE IS MUTEXED!!! */

// pthread_mutex_unlock( &mutex );
// /* END MUTEX HERE     **/


/*******************************/
 
 
double g_time_in_secs = 0;
double g_processor_frequency = 1600000000.0; // processing speed for BG/Q
unsigned long long g_start_cycles=0;
unsigned long long g_end_cycles=0;

// You define these

struct rank_data{
    short*** board;
    short** ghost_above;
    short** ghost_below;
    int* i;
};
int rows_per_thread;
int rows_per_rank;
int my_rank;
int num_ranks;
int num_threads = 2;

/***************************************************************************/
/* Function Decs ***********************************************************/
/***************************************************************************/

// You define these
short** make_board(int rows);
void deallocate_mem(short*** arr, int rows);
void print_board(short** board, int rows);
short* make_ghost_row();
short** copy_board(short** board, int rows);
void exchange_ghosts(struct rank_data * x);
void* thread_init(void*);
short** copy_board_with_start(short** board, int rows, int start);
void print_row(short* row);
void copy_thread_rows_to_universe(short*** board, short*** thread_board, int rows, int start);
void copy_ghost_row_to_universe(short** ghost, short** thread_ghost);
short* copy_row(short* ghost_row);

/***************************************************************************/
/* Function: Main **********************************************************/
/***************************************************************************/

int main(int argc, char *argv[])
{
// Example MPI startup and using CLCG4 RNG
    MPI_Init( &argc, &argv);
    MPI_Comm_size( MPI_COMM_WORLD, &num_ranks);
    MPI_Comm_rank( MPI_COMM_WORLD, &my_rank);
    
// Init 32,768 RNG streams - each rank has an independent stream
    InitDefault();
    
// Note, used the my_rank to select which RNG stream to use.
// You must replace my_rank with the right row being used.
// This just show you how to call the RNG.    
    // printf("Rank %d of %d has been started and a first Random Value of %lf\n", 
	   // my_rank, num_ranks, GenVal(my_rank));
    
    MPI_Barrier( MPI_COMM_WORLD );
    
// Insert your code
    // Following algorithm description from class
    // =========================================================
    // =========================================================
    // =========================================================
    // =========================================================
    // =========================================================
    

    // if rank 0 / pthread0, start time with GetTimeBase() 
    if(my_rank == 0){
        g_start_cycles = GetTimeBase();
        
        // sanity check
        printf("Board size: %d x %d\n", board_size, board_size);
        int cells_per_rank = board_size*(board_size/num_ranks);
        int rows_per_rank = board_size/num_ranks;
        printf("each of %d ranks will have %d cells, %d rows, %d ghost rows\n",
            num_ranks,cells_per_rank, rows_per_rank, 2);
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
    int num_ticks = 1;
    rows_per_rank = board_size/num_ranks; 
    rows_per_thread = rows_per_rank/num_threads;

    short** board = make_board(rows_per_rank);
    short* ghost_above = make_ghost_row();
    short* ghost_below = make_ghost_row(); 

    if(my_rank==0){
        board[0][1] = 0;
        board[rows_per_rank-1][2] = 2;
    }

    if(my_rank==num_ranks-1){
        board[0][0] = 3;
        board[rows_per_rank-1][board_size-1] = 4;
    }
  
    /*  3
        For all number of ticks, complete a round of the GOL
    */
    struct rank_data thread_data;
    for(int tick=0; tick<num_ticks; tick++){
        //print_board(board);
        pthread_t tid[num_threads];
        for(int i=0;i<num_threads;i++){
            int thread_val = i;

            thread_data.ghost_above = &ghost_above;
            thread_data.ghost_below = &ghost_below;
            thread_data.board = &board;
            thread_data.i = malloc(sizeof(int *));
            *(thread_data.i) = i;

            // MPI_Barrier(MPI_COMM_WORLD);
            // for(int i = 0; i < num_ranks; ++i){
            //     if(i== my_rank){
            //         printf("Rank: %d board\n", i);
            //         print_board(*(thread_data.board), rows_per_thread);
            //         printf("\n");
            //     }
            // }
            // MPI_Barrier(MPI_COMM_WORLD);
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
            // first thread so send and receive ghost row above
            
            // if(my_rank==num_ranks-1 && i==num_threads-1){
            //     printf("\nlast rank ghost row below:\n");
            //     print_row(*(thread_data.ghost_below));
            //     printf("\nlast rank ghost row above:\n");
            //     print_row(*(thread_data.ghost_above));
            // }

            
            /*  
                Create Pthreads here. All threads should go into for loop
            */
            int rc = pthread_create(&tid[i], NULL, thread_init, &thread_data);
            if (rc != 0) {
                fprintf(stderr, "ERROR: pthread_create() failed\n");
            }

            pthread_join(tid[i], NULL);
        }

        printf("\n Rank %d ghost row below:\n", my_rank);
        print_row(ghost_below);
        printf("\n Rank %d ghost row above:\n", my_rank);
        print_row(ghost_above);

        


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
    
    // for (int i = 0; i<rows_per_thread; i++){
    //     free(board[i]);
    // }
    // free(board); 
    // free(ghost_above);
    // free(ghost_below);

    
    

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
short** make_board(int rows){
    short** board = calloc(rows,sizeof(short*));
    for(int i=0;i<rows;++i){
        board[i] = calloc(board_size,sizeof(short));
    }

    for(int i=0;i<rows;++i){
        for(int j=0;j<board_size;++j){
            board[i][j] = ALIVE;
        }
    }
    return board;   
}

short* make_ghost_row(){
    short* ghost = calloc(board_size,sizeof(short));
    for(int i=0;i<board_size;++i){
        ghost[i] = ALIVE;
    }
    return ghost;
}

short** copy_board(short** board, int rows){
    short** new_board = make_board(rows);
    for(int i=0;i<rows;++i){
        for(int j=0;j<board_size;++j){
            new_board[i][j] = board[i][j];
        }
    }
    return new_board;
}

short* copy_row(short* ghost_row){
    short* new_ghost_row = make_ghost_row();
    for(int i=0;i<board_size;++i){
        new_ghost_row[i] = ghost_row[i];
        
    }
    return new_ghost_row;
}

short** copy_board_with_start(short** board, int rows, int start){
    short** new_board = make_board(rows);
    for(int i=0;i<rows;++i){
        for(int j=0;j<board_size;++j){
            new_board[i][j] = board[i+start][j];
        }
    }
    return new_board;
}

void copy_thread_rows_to_universe(short*** board, short*** thread_board, int rows, int start){
    for(int i=0;i<rows;++i){
        for(int j=0;j<board_size;++j){
            *board[i+start][j] = *thread_board[i][j];
        }
    }
}

void copy_ghost_row_to_universe(short** ghost, short** thread_ghost){
    for(int i=0;i<board_size;++i){
        *ghost[i] = *thread_ghost[i];   
    }
}

void deallocate_mem(short*** arr, int rows){
    for (short i = 0; i<rows; i++){
        free((*arr)[i]);
    }
    free(*arr); 
}

void print_board(short** board, int rows){
    for(int i=0;i<rows;++i){
        for(int j=0;j<board_size;++j){
            printf("%hd",board[i][j]);
        }
        printf("\n");      
    }
}

void print_row(short* row){
    for(int i=0;i<board_size;++i){
        printf("%hd",row[i]); 
    }
    printf("\n"); 
}

void exchange_ghosts(struct rank_data * x){
    struct rank_data thread_data = *x;

    MPI_Barrier(MPI_COMM_WORLD);
    if(num_ranks>1 && *(thread_data.i)==0){
        // tag of 0
        MPI_Request send_request, recv_request;
        MPI_Status status;
        
        if(my_rank==0){
            //very top ghost above receives from very last row 
            MPI_Irecv(*(thread_data.ghost_below), board_size, MPI_SHORT, num_ranks-1, 0, MPI_COMM_WORLD, &recv_request);
            //very top row sends to ghost below of last rank
            MPI_Isend((*(thread_data.board))[0], board_size, MPI_SHORT, num_ranks-1, 1, MPI_COMM_WORLD, &send_request);
            
            printf("before 1\n");
            printf("%d receiving from %d\n", my_rank, num_ranks-1);
            MPI_Wait(&recv_request, &status);
            printf("after 1\n");
            printf("before 2\n"); 
            printf("%d sending to %d\n", my_rank, num_ranks-1);
            MPI_Wait(&send_request, &status);    
            printf("after 2\n");  

            // printf("Rank %d, thread %d, receiving ghost row below from rank %d\n",my_rank, *(thread_data.i), num_ranks-1);
            // print_row(*(thread_data.ghost_below));
        }

        else{
            //every normal ghost above receives from rank-1's last row
            MPI_Irecv(*(thread_data.ghost_below), board_size, MPI_SHORT, my_rank-1, 1, MPI_COMM_WORLD, &recv_request);
            //every top row above sends to rank-1's ghost below
            MPI_Isend((*(thread_data.board))[0], board_size, MPI_SHORT, my_rank-1, 0, MPI_COMM_WORLD, &send_request);
            printf("before 3\n");
            printf("%d receiving from %d\n", my_rank, my_rank-1);
            MPI_Wait(&recv_request, &status); 
            printf("after 3\n");
            printf("before 4\n");
            printf("%d sending to %d\n", my_rank, my_rank-1);
            MPI_Wait(&send_request, &status);
            printf("after 4\n");

            // printf("Rank %d, thread %d, receiving ghost row below from rank %d\n",my_rank, *(thread_data.i), my_rank-1);
            // print_row(*(thread_data.ghost_below));
        }
    }
    // last thread so send and receive ghost rows
    if(num_ranks>1 && *(thread_data.i)==num_threads-1){
        // tag of 1
        MPI_Request send_request, recv_request;
        MPI_Status status;
        
        if(my_rank==num_ranks-1){
            //very bottom ghost below receives from very top row
            MPI_Irecv(*(thread_data.ghost_above), board_size, MPI_SHORT, 0, 1, MPI_COMM_WORLD, &recv_request);
            //very bottom row sends to ghost above of first rank
            MPI_Isend((*(thread_data.board))[rows_per_rank-1], board_size, MPI_SHORT, 0, 0, MPI_COMM_WORLD, &send_request);
            
            printf("before 5\n");
            printf("%d receiving from %d\n", my_rank, 0);
            MPI_Wait(&recv_request, &status); 
            printf("after 5\n");
            printf("before 6\n");
            printf("%d sending to %d\n", my_rank, 0);
            MPI_Wait(&send_request, &status);
            printf("after 6\n");

            // printf("Rank %d, thread %d, receiving ghost row above from rank %d\n",my_rank, *(thread_data.i), 0);
            // print_row(*(thread_data.ghost_above));

        }

        else{
            //every normal ghost below receives form rank+1's first row
            MPI_Irecv(*(thread_data.ghost_above), board_size, MPI_SHORT, my_rank+1, 0, MPI_COMM_WORLD, &recv_request);
            //every normal bottom row sends to rank+1's ghost above
            MPI_Isend((*(thread_data.board))[rows_per_rank-1], board_size, MPI_SHORT, my_rank+1, 1, MPI_COMM_WORLD, &send_request);
            
            printf("before 7\n");
            printf("%d receiving from %d\n", my_rank, my_rank+1);
            MPI_Wait(&recv_request, &status); 
            printf("after 7\n");
            printf("before 8\n");
            printf("%d sending to %d\n", my_rank, my_rank+1);
            MPI_Wait(&send_request, &status);
            printf("after 8\n");

            // printf("Rank %d, thread %d, receiving ghost row above from rank %d\n",my_rank, *(thread_data.i), my_rank+1);
            // print_row(*(thread_data.ghost_above));
        }
    }
}

void* thread_init(void* x){
    struct rank_data thread_data = *((struct rank_data *)x);

    exchange_ghosts(&thread_data);

    pthread_exit(NULL);
}











