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
 
___________________________________________________________

4/8/19 
Issues - deadlock when used with more than 2 mpi ranks
            - set DEBUG to 1 to try and fix that issue
            - DEBUG is only used with this issue
       - write to file doesnt actually write anything to the file
            - just get nothing

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

struct thread_struct{
    short*** board;
    short** ghost_above;
    short** ghost_below;
    int* i;
    int* current_tick;
};
int rows_per_thread;
int rows_per_rank;
int my_rank;
int num_ranks;
int num_threads = 8;
int DEBUG = 1;
int PRINT = 0;
float threshold = 0.5;
#define board_size 64//32768
int num_ticks = 1;
long long * tick_sums;
int num_bins = 8;


/***************************************************************************/
/* Function Decs ***********************************************************/
/***************************************************************************/

// You define these
short** make_board(int rows);
void deallocate_mem(short*** arr, int rows);
void print_board(short** board, int rows);
short* make_ghost_row();
short** copy_board(short** board, int rows);
void exchange_ghosts(struct thread_struct * x);
void* thread_init(void*);
short** copy_board_with_start(short** board, int rows, int start);
void print_row(short* row);
void copy_thread_rows_to_universe(short*** board, short*** thread_board, int rows, int start);
void copy_ghost_row_to_universe(short** ghost, short** thread_ghost);
short* copy_row(short* ghost_row);
void write_to_file(short** board, int h, int w);
void bin_board(struct thread_struct * x);

/***************************************************************************/
/* Function: Main **********************************************************/
/***************************************************************************/

int main(int argc, char *argv[])
{
// Example MPI startup and using CLCG4 RNG
    MPI_Init( &argc, &argv);
    MPI_Comm_size( MPI_COMM_WORLD, &num_ranks);
    MPI_Comm_rank( MPI_COMM_WORLD, &my_rank);
    
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

        tick_sums = calloc(num_ticks,sizeof(long long));
    }
    
    rows_per_rank = board_size/num_ranks; 
    rows_per_thread = rows_per_rank/num_threads;

    short** board = make_board(rows_per_rank);
    short* ghost_above = make_ghost_row();
    short* ghost_below = make_ghost_row(); 

    if(num_ranks==1){
        for(int i=0;i<board_size;i++){
            ghost_above[i] = board[board_size-1][i];
            ghost_below[i] = board[0][i];
        }
    }
  
    struct thread_struct thread_data;
    for(int tick=0; tick<num_ticks; tick++){
        if(PRINT){
            printf("Generation %d\n",tick+1);    
        }
        
        pthread_t tid[num_threads];
        for(int i=0;i<num_threads;i++){
            InitDefault();
            thread_data.ghost_above = &ghost_above;
            thread_data.ghost_below = &ghost_below;
            thread_data.board = &board;
            thread_data.i = malloc(sizeof(int *));
            *(thread_data.i) = i;
            thread_data.current_tick = malloc(sizeof(int *));
            *(thread_data.current_tick) = tick;
          
            int rc = pthread_create(&tid[i], NULL, thread_init, &thread_data);
            if (rc != 0) {
                fprintf(stderr, "ERROR: pthread_create() failed\n");
            }
            pthread_join(tid[i], NULL);
        }    
    }
    g_end_cycles = GetTimeBase();

    MPI_Barrier(MPI_COMM_WORLD);
    write_to_file(board, board_size, board_size);

    //bin_board(&thread_data);

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
short** make_board_cols(int rows, int cols){
    short** board = calloc(rows,sizeof(short*));
    for(int i=0;i<rows;++i){
        board[i] = calloc(cols,sizeof(short));
    }

    for(int i=0;i<rows;++i){
        for(int j=0;j<cols;++j){
            board[i][j] = 0;
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
    MPI_Barrier(MPI_COMM_WORLD);
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

void print(struct thread_struct * x){
    struct thread_struct thread_data = *x;
    print_board((*(thread_data.board)),rows_per_rank);
    //printf("\n");
}

void write_to_file(short** board, int h, int w){
    // if(my_rank==0){
    //     MPI_Status status;
    //     MPI_File file;
    //     MPI_File_open(MPI_COMM_SELF, "test.txt", MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &file);
    //     MPI_File_write(file,*board,4,MPI_SHORT,&status);
    //     MPI_File_close(&file);
    // }
    // return;
    printf("Rank %d writing to file...\n",my_rank);
    h = h/num_ranks;

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Status status;
    MPI_File file;
    MPI_File_open(MPI_COMM_SELF, "./test.txt", MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &file);
    
    // writing each row individually because the write buffer wants a void*
    MPI_Barrier(MPI_COMM_WORLD);
    printf("w is %d\n",w);
    for(int i=0;i<h;i++){

        // printf("Rank %d thread %d writing to file at byte offset\n",my_rank, i);
        // for(int j=0;j<w-1;j+=2){
        //     MPI_Offset offset = ((my_rank * h + i) * 2*w +j) * sizeof(short);
        //     char c[2];
        //     c[0] = (board[i][j]+'0');
        //     c[1] = ',';
        //     MPI_File_write_at(file, offset, c, 2*sizeof(short), MPI_SHORT, &status);
        // }

        // ==================================
        int line_alloc = w+2;
        char string[line_alloc];
        for(int j=0;j<w;j++){
            string[j] = 1+'0';
        }
        string[line_alloc - 2] = '\n';
        string[line_alloc - 1] = '\0';
        MPI_Offset offset = (my_rank * h  + i)*line_alloc;
        printf("Rank %d offset %lld printing %s\n",my_rank,offset,string);
        MPI_File_write_at(file, offset, string, strlen(string), MPI_CHAR, &status);
        //MPI_File_write(file,string,strlen(string),MPI_CHAR,&status);
    }    

    MPI_File_close(&file);
    printf("Finished writing to file.\n");
}

void exchange_ghosts(struct thread_struct * x){
    struct thread_struct thread_data = *x;
    MPI_Request send_request, recv_request;
    MPI_Status status;
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
    MPI_Barrier(MPI_COMM_WORLD);

    if(*(thread_data.i)==0){
        if(my_rank>0){
            MPI_Irecv(*(thread_data.ghost_above), board_size, MPI_SHORT, my_rank-1, 1, MPI_COMM_WORLD, &recv_request);
            MPI_Wait(&recv_request, &status); 
        }
        MPI_Isend((*(thread_data.board))[rows_per_rank-1], board_size, MPI_SHORT, (my_rank+1)%num_ranks, 1, MPI_COMM_WORLD, &send_request);
        MPI_Wait(&send_request, &status);
        if(my_rank==0){
            MPI_Irecv(*(thread_data.ghost_above), board_size, MPI_SHORT, num_ranks-1, 1, MPI_COMM_WORLD, &recv_request);
            MPI_Wait(&recv_request, &status); 
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    if(*(thread_data.i)==num_threads-1){
        if(my_rank<num_ranks-1){
            MPI_Irecv(*(thread_data.ghost_below), board_size, MPI_SHORT, my_rank+1, 0, MPI_COMM_WORLD, &recv_request);
            MPI_Wait(&recv_request, &status);
        }
        int send_to = (my_rank-1)%num_ranks;
        if(send_to==-1){
            send_to = num_ranks-1;
        }
        MPI_Isend((*(thread_data.board))[0], board_size, MPI_SHORT, send_to, 0, MPI_COMM_WORLD, &send_request);
        MPI_Wait(&send_request, &status);
        if(my_rank==num_ranks-1){
            MPI_Irecv(*(thread_data.ghost_below), board_size, MPI_SHORT, 0, 0, MPI_COMM_WORLD, &recv_request);
            MPI_Wait(&recv_request, &status); 
        }
    }
    // ==================================================
}

int get_state(struct thread_struct * x, int row, int col){
    struct thread_struct thread_data = *x;
    if(col==-1){
        return get_state(x,row,board_size - 1);
    }
    if(col==board_size){
        return get_state(x,row,0);
    }
    
    if(row==-1){
        //printf("(%d,%d) %hd\n",row,col,(*(thread_data.ghost_above))[col]);
        return (*(thread_data.ghost_above))[col];
    }
    if(row==rows_per_rank){
        //printf("(%d,%d) %hd\n",row,col,(*(thread_data.ghost_below))[col]);
        return (*(thread_data.ghost_below))[col];
    }
    //printf("(%d,%d) %hd\n",row,col,(*(thread_data.board))[row][col]);
    return (*(thread_data.board))[row][col];
}

int find_num_neighbors(struct thread_struct * x, int row, int col){
    int num_neighbors=0;
    // row above
        // top left
        num_neighbors += get_state(x,row-1,col-1);
        // top 
        num_neighbors += get_state(x,row-1,col);
        // top right
        num_neighbors += get_state(x,row-1,col+1);
    // same row
        // left
        num_neighbors += get_state(x,row,col-1);
        // right
        num_neighbors += get_state(x,row,col+1);
    // row below
        // bottom left
        num_neighbors += get_state(x,row+1,col-1);
        // bottom
        num_neighbors += get_state(x,row+1,col);
        // bottom right
        num_neighbors += get_state(x,row+1,col+1);
    //printf("(%d,%d) diag %d\n",row,col,get_state(x,row-1,col-1));
    return num_neighbors;
}

void apply_rules(struct thread_struct * x){
    struct thread_struct thread_data = *x;
    /*

    Rules:
        Any live call with fewer than two live neighbors dies
        Any live cell with more than 3 live neighbors dies
        Any other live cell lives on

        Any dead cell with exactly 3 neighbors lives

    */

    for(int i=0; i<rows_per_rank; i++){
        int global_index = i+(*(thread_data.i)*rows_per_rank);
        for(int j=0; j<board_size; j++){
            if(PRINT==2){
                print(x);
            }
            float val = GenVal(global_index);
            //printf("rank %d row %d: %f\n", my_rank, global_index, val);
            
            // apply rules
            if(val > threshold){

                int num_neighbors = find_num_neighbors(x,i,j);
                if(PRINT==2){
                    printf("(%d,%d) neighbors: %d\n", i, j, num_neighbors);
                }
                // Any live cell with fewer than two live neighbors dies
                if(num_neighbors<2){
                    (*(thread_data.board))[i][j] = DEAD;
                }
                // Any live cell with more than 3 live neighbors dies
                if(num_neighbors>3){
                    (*(thread_data.board))[i][j] = DEAD;
                }
                // Any dead cell with exactly 3 neighbors lives
                if(num_neighbors==3){
                    (*(thread_data.board))[i][j] = ALIVE;
                }
            }

            // pick random state of LIVE or DEAD
            else{
                float state_rand = GenVal(global_index);
                if(PRINT==2){
                    printf("(%d,%d) Random state: %f\n",i,j,state_rand);
                }
                
                if(state_rand>0.5){
                    (*(thread_data.board))[i][j] = DEAD;
                }
                else{
                    (*(thread_data.board))[i][j] = ALIVE;   
                }
            }
        }
    }
}

long long sum_rank(struct thread_struct * x){
    struct thread_struct thread_data = *x;
    long long total = 0;
    for(int i=0; i<rows_per_rank; i++){
        for(int j=0; j<board_size; j++){
            total += (*(thread_data.board))[i][j];
        }
    }
    return total;
}

void bin_board(struct thread_struct * x){
    struct thread_struct thread_data = *x;
    printf("%d\n",*(thread_data.i));

    // the number of cells being summed into each bin
    int length_of_bin = board_size/num_bins;

    // the vertical number of bins going in each rank
    int num_vertical_bins = num_bins/num_threads;

    // for each row in each rank
    for(int row=0;row<rows_per_rank;row++){
        for(int bin=0;bin<num_bins;bin++){

        }
    }

}

void* thread_init(void* x){
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

    struct thread_struct thread_data = *((struct thread_struct *)x);

    exchange_ghosts(&thread_data);
    if(PRINT){
        MPI_Barrier(MPI_COMM_WORLD);
        printf("\n");
        print(x);
    }

    // test ghost row exchange
    // MPI_Barrier(MPI_COMM_WORLD);
    // if(my_rank==0){
    //     printf("\n");
    //     print_row(*(thread_data.ghost_above));
    //     print_row(*(thread_data.ghost_below));
    // }

    apply_rules(&thread_data);
    if(PRINT){
        MPI_Barrier(MPI_COMM_WORLD);
        printf("\n");
        print(x);
    }
    
    long long rank_sum = sum_rank(x);
    MPI_Reduce(&rank_sum, &tick_sums[*(thread_data.current_tick)], 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    if(my_rank==0 && PRINT){
        printf("%lld alive on board\n\n", tick_sums[*(thread_data.current_tick)]);
    }

    pthread_exit(NULL);
    return NULL;
}











