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

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// CONTROL AREA
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
int num_threads = 1;
int DEBUG = 0;
int PRINT = 0;
float threshold = 0.25;
#define board_size 32768
int num_ticks = 256;
long long * tick_sums;
int num_bins = 1024;
int io = 1;
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// DONT EDIT BEYOND HERE
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

// for output on BGmfQ
char* master_filename_title = "./Config_1.txt";
FILE * master_file;



/***************************************************************************/
/* Function Decs ***********************************************************/
/***************************************************************************/

// You define these
short** make_board(int rows);
void deallocate_mem(short*** arr, int rows);
void print_board(short** board, int rows);
short* make_ghost_row();
void exchange_ghosts(struct thread_struct * x);
void* thread_init(void*);
void print_row(short* row);
void write_to_file(short** board, int h, int w);

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
    
    // master thread does no computation, its all done in thread init
    if(num_threads==0){
        num_threads+=1;
    }

    // if rank 0 / pthread0, start time with GetTimeBase() 
    if(my_rank == 0){
        g_start_cycles = GetTimeBase();
        tick_sums = calloc(num_ticks,sizeof(long long));
        printf("Opening %s\n",master_filename_title);
        master_file = fopen(master_filename_title, "w");
        fprintf(master_file,"%s","~~~Begin File~~~\n");
        printf("Should have written\n");
        fprintf(master_file, "\nBoard is %dx%d\n", board_size,board_size);
        fprintf(master_file, "Number of Ranks: %d\n", num_ranks);
        fprintf(master_file, "Number of Threads: %d\n", num_threads);
        fprintf(master_file, "Number of Ticks: %d\n", num_ticks);
        fprintf(master_file, "Threshold set to %.2f\n", threshold);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    
    // set some helpful global vars
    rows_per_rank = board_size/num_ranks; 
    rows_per_thread = rows_per_rank/num_threads;

    // initialize the board for each rank, these will be 
    // the only board and ghost allocations
    short** board = make_board(rows_per_rank);
    short* ghost_above = make_ghost_row();
    short* ghost_below = make_ghost_row(); 

    // if only one master rank, ghost rows are needed still 
    // but are never reset after this point
    if(num_ranks==1){
        for(int i=0;i<board_size;i++){
            ghost_above[i] = board[board_size-1][i];
            ghost_below[i] = board[0][i];
        }
    }

    struct thread_struct thread_data;
    for(int tick=0; tick<num_ticks; tick++){
        if(my_rank==0){
            printf("tick: %d\n",tick);
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

    if(io){
		if(my_rank==0) printf("Begin Parallel Write\n");
        unsigned long long start_io_time = GetTimeBase();
        write_to_file(board, board_size, board_size);
        unsigned long long end_io_time = GetTimeBase();
	
        if(my_rank==0){
        	printf("End Parallel Write\n");
            int time_in_secs = ((double)(end_io_time - start_io_time));// / g_processor_frequency;
            fprintf(master_file, "\nParallel I/O time: %ds\n\n", time_in_secs);      
        }
        
    }
    else{
        if(my_rank==0){
            fprintf(master_file, "\n%s\n\n", "I/O is turned off for this configuration");
        }
    }

    if(my_rank==0){
        fprintf(master_file,"Alive after each generation:\n["); 
        for(int i=0;i<num_ticks;i++){
            fprintf(master_file,"%lld, ", tick_sums[i]);
        } 
        fprintf(master_file,"]\n\n");
        
        int time_in_secs = ((double)(g_end_cycles - g_start_cycles));// / g_processor_frequency;
        fprintf(master_file,"Total execution time: %ds\n",time_in_secs);
        fclose(master_file);
    }

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
}

void write_to_file(short** board, int h, int w){
    h = h/num_ranks;
    w = w+2;

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Status status;
    MPI_File file;
    MPI_File_open(MPI_COMM_SELF, "./test_output.txt", MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &file);
    
    // writing each row individually because the write buffer wants a void*
    MPI_Barrier(MPI_COMM_WORLD);
    for(int j=0;j<h;j++){
        char string[w];
        for(int i=0;i<w-1;i++){
            string[i] = board[j][i]+'0';
        }
        string[w-1] = '\n';
        MPI_Offset offset = (my_rank * h  + j)*(w);
        //printf("Rank %d offset %lld printing %s\n",my_rank,offset,string);
        MPI_File_write_at(file, offset, string, w, MPI_CHAR, &status);
    }    
    MPI_File_close(&file);
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
    

    pthread_exit(NULL);
    return NULL;
}











