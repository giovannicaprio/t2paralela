/*

Para compilar no MAC 
gcc-8 -fopenmp paralelo.c -o paralelo

mpicc mpi_paralelo.c -o mpi_paralelo
mpirun -np 2 ./mpi_paralelo

*/

#include "mpi.h"

#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <sys/time.h>
#define k 40
#define position(row,col) //printf("%c[%d;%dH",k,row,col)
#define clear() //printf("%c[2J",k)
#define ROWS 175

char  *TOP[ROWS];
char  *BOTTOM[ROWS];
char  *FULL[ROWS] = {
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "....O...................................",
  ".....O..................................",
  "...OOO..................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "........................................",
  "..........O.............................",
  "...........O............................",
  ".........OOO............................",
};

int ROWSIZE = strlen("........................................") + 1;

/***
GLOBAL MPI variables 
***/

int my_rank;  /* Identificador do processo */
int proc_n;   /* Número de processos */
int source;   /* Identificador do proc.origem */
int dest;     /* Identificador do proc. destino */
int tag = 50; /* Tag para as mensagens */
   
char message[100]; /* Buffer para as mensagens */
MPI_Status status; /* Status de retorno */


void game( char**, char**);
void gameMPI( char**, char**);
void initialize( );
void display( char ** );
void life( char**, char**, int );



// --------------------------------------------------------------------------
// initialize
// initialize top and bottom of grid
void initialize( ) {
  int x;
  //allocate memory and copy input grid to top and bottom for processing
  for (x = 0; x< ROWS; x++ )  {
    TOP[x] = (char *) malloc( (strlen( FULL[0] ) + 1 ) * sizeof( char ) );
    strcpy( TOP[x], FULL[x] );
    BOTTOM[x] = (char *) malloc( (strlen( TOP[0] )+1) * sizeof( char )  );
    strcpy( BOTTOM[x], TOP[x] );
  }
}

// Display current iteration
void display( char* dish[] ) {
  int j;
  for (j=0; j<ROWS; j++ ) {
    position( j, 0 );
    printf( "%s \n", dish[j] );
  }

}

void game( char** part, char** nextStep) {
  /*
   * When given a portion of the grid, this function will calculate the 
   * next step and return the new portion of the grid
   */
  int m, n, row;
  int rowSize = strlen( part[0] );
  int partSize = ROWS;

  int th_id, nthreads;

  #pragma omp parallel for private(row, m) schedule(static)
  for (row = 0; row < ROWS; row++) {            // loop through the rows in grid
    for ( m = 0; m < rowSize; m++) {            // loop through each character in the row
      int o, n, side = 0;
      char curr = part[row][m];
        // search a 3 by 3 area in the grid for any live cells 'O'
        for ( o = row - 1; o <= row + 1; o++) {
           int realo = o;// traverse from bottom to top
           if (o == -1)
             realo = partSize - 1;
           if (o == partSize)
             realo = 0;
           for ( n = m - 1; n <= m + 1; n++) {
             int realn = n; // traverse left to right
             if (n == -1)
               realn = rowSize - 1;
             if (n == rowSize)
               realn = 0;
             if (o == row && n == m)
               continue;
             if (part[realo][realn] == 'O')
               side++;
           }
        }
      //implementacao da regra do jogo de acordo com o que encontrou em volta
      if (curr == 'O') {
        if (side < 2 || side > 3)
          nextStep[row][m] = '.';
        else
          nextStep[row][m] = 'O';
      }
      if (curr == '.') {
        if (side == 3)
          nextStep[row][m] = 'O';
        else
          nextStep[row][m] = '.';
      }
    }
  }
}

void gameMPI( char** part, char** nextStep) {
  /*
   * When given a portion of the grid, this function will calculate the 
   * next step and return the new portion of the grid
   */
  int m, n, row;
  int rowSize = strlen( part[0] );
  int partSize = ROWS;

  int th_id, nthreads;

    /*
    if(my_rank == 0){ //master
       //printf("rank %d out of %d processors\n",my_rank, proc_n);
      MPI_Send ("1", strlen(message)+1, MPI_CHAR,1, tag, MPI_COMM_WORLD);
      
      MPI_Recv (message, 100, MPI_CHAR , 1, tag, MPI_COMM_WORLD, &status);
      printf("%s\n", message);


    }else{ //slaves
      MPI_Recv (message, 100, MPI_CHAR , source, tag, MPI_COMM_WORLD, &status);
      printf("%s\n", message);
      MPI_Send ("2", strlen(message)+1, MPI_CHAR,0, tag, MPI_COMM_WORLD);

    }
    */

  #pragma omp parallel for private(row, m) schedule(static)
  for (row = 0; row < ROWS; row++) {            // loop through the rows in grid
    for ( m = 0; m < rowSize; m++) {            // loop through each character in the row
      if(my_rank == 0){ //master
        //dispara os escravos passando 
        MPI_Send (&row, strlen(message)+1, MPI_CHAR,1, tag, MPI_COMM_WORLD);
    
        //recebe a reposta dos escravos
        MPI_Recv (message, 100, MPI_CHAR , 1, tag, MPI_COMM_WORLD, &status);
        printf("%s\n", message);
      }else{
        //recebe a tarefa do master
         MPI_Recv (message, 100, MPI_CHAR , source, tag, MPI_COMM_WORLD, &status);
         printf("row eh%s\n", message);


        //devolve o resultado pro master
        MPI_Send ("2", strlen(message)+1, MPI_CHAR,0, tag, MPI_COMM_WORLD);



      }



      int o, n, side = 0;
      char curr = part[row][m];
        
        // search a 3 by 3 area in the grid for any live cells 'O'
        for ( o = row - 1; o <= row + 1; o++) {
           int realo = o;// traverse from bottom to top
           if (o == -1)
             realo = partSize - 1;
           if (o == partSize)
             realo = 0;
           for ( n = m - 1; n <= m + 1; n++) {
             int realn = n; // traverse left to right
             if (n == -1)
               realn = rowSize - 1;
             if (n == rowSize)
               realn = 0;
             if (o == row && n == m)
               continue;
             if (part[realo][realn] == 'O')
               side++;
           }
        }
      //implementacao da regra do jogo de acordo com o que encontrou em volta
      if (curr == 'O') {
        if (side < 2 || side > 3)
          nextStep[row][m] =  '.';
        else
          nextStep[row][m] = 'O';
      }
      if (curr == '.') {
        if (side == 3)
          nextStep[row][m] = 'O';
        else
          nextStep[row][m] = '.';
      }
    }
  }
}


//function to get wall time
double get_wall_time(){
  struct timeval time;
  if (gettimeofday(&time,NULL)){
    //  Handle error
    return 0;
  }
  return (double)time.tv_sec + (double)time.tv_usec * .000001;
}

// --------------------------------------------------------------------------
void life( char** part, char** nextPart, int my_rank ) {

  /*
   * Given an array of string representing the current population of cells
   * in a petri dish, computes the new generation of cells according to
   * the rules of the game. A new array of strings is returned.
   */
  int i, j, row;
  int rowLength = strlen(part[0]);
  //printf("\r rowLength rowLength orowLength: %d\n\n", rowLength); 
  //printf("\r quel eh o gaucho que ta aqui: %d\n\n", my_rank); 

  int dishLength = ROWS;
  
  int lowerRow, upperRow;

  //--- slice the array into two halves.  Rank 0 is lower half, ---
  //--  Rank 1 is upper half.                                   ---
  if ( my_rank == 0 ) {
    lowerRow = 0;
    upperRow = ROWS/2;
  }
  if ( my_rank == 1 ) {

    lowerRow = ROWS/2;
    upperRow = ROWS;
  }
  
  for (row = lowerRow; row < upperRow; row++) {// each row
    
    if (part[row] == NULL )
      continue;

    for ( i = 0; i < rowLength; i++) { // each char in the
                                       // row

      int r, j, neighbors = 0;
      char current = part[row][i];

      // loop in a block that is 3x3 around the current cell
      // and count the number of 'O' cells.
      for ( r = row - 1; r <= row + 1; r++) {

        // make sure we wrap around from bottom to top
        int realr = r;
        if (r == -1)
          realr = dishLength - 1;
        if (r == dishLength)
          realr = 0;

        for (int j = i - 1; j <= i + 1; j++) {

          // make sure we wrap around from left to right
          int realj = j;
          if (j == -1)
            realj = rowLength - 1;
          if (j == rowLength)
            realj = 0;

          if (r == row && j == i)
            continue; // current cell is not its
          // neighbor
          if (part[realr][realj] == 'O')
            neighbors++;
        }
      }

      if (current == 'O') {
        if (neighbors < 2 || neighbors > 3)
          nextPart[row][i] = '.';
        else
          nextPart[row][i] = 'O';
      }

      if (current == '.') {
        if (neighbors == 3)
          nextPart[row][i] = 'O';
        else
          nextPart[row][i] = '.';
      }
    }
  }
}



int main( int argc, char* argv[] ) {
  /************
      MPI
  ************/
  double t1,t2;
  t1 = MPI_Wtime();  // inicia a contagem do tempo  


  double start,end;
  start = get_wall_time();
  int steps = 10;      // # of steps in game
  int u;
  char **part, **nextPart, **tmp;
  //empty screen for next step
  clear();
  initialize();
  part   = TOP;
  nextPart = BOTTOM;

  MPI_Init (&argc , & argv);

  char processor_name[MPI_MAX_PROCESSOR_NAME];

  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &proc_n);

  for (u = 0; u < steps; u++) { // show rest of the steps
    // follow rules in game of life to current step and display next step
    //game(part, nextPart);
    //gameMPI(part, nextPart);

    printf("rank %d out of %d proc_n\n",my_rank, proc_n);

    life(part, nextPart, my_rank);

    //display( part );
    //printf("%s\n"," --- NOVA JOGADA ---" );
    // This is just for user to be able to see whats going on. This must be removed for timing calculations
    //sleep(1); 
    

    if (my_rank == 0 ) {
      //          buffer                #items   item-size src/dest tag   world  
      MPI_Send( nextPart[    0   ], ROWSIZE, MPI_CHAR,    1,     0,  MPI_COMM_WORLD );
      MPI_Send( nextPart[ROWS/2-1], ROWSIZE, MPI_CHAR,    1,     0,  MPI_COMM_WORLD );
      MPI_Recv( nextPart[ROWS-1],   ROWSIZE, MPI_CHAR,    1,     0,  MPI_COMM_WORLD, &status );
      MPI_Recv( nextPart[ROWS/2],   ROWSIZE, MPI_CHAR,    1,     0,  MPI_COMM_WORLD, &status );
    }
    if (my_rank == 1 ) {
      MPI_Recv( nextPart[ 0      ], ROWSIZE, MPI_CHAR,    0,     0,  MPI_COMM_WORLD, &status );
      MPI_Recv( nextPart[ROWS/2-1], ROWSIZE, MPI_CHAR,    0,     0,  MPI_COMM_WORLD, &status );
      MPI_Send( nextPart[ROWS-1],   ROWSIZE, MPI_CHAR,    0,     0,  MPI_COMM_WORLD );
      MPI_Send( nextPart[ROWS/2],   ROWSIZE, MPI_CHAR,    0,     0,  MPI_COMM_WORLD );
    }


    tmp = part; // obtain next step
    part = nextPart;
    nextPart = tmp;
  }

  end = get_wall_time();
  //printf("Total time elasped is %lf \n",(end-start));
  t2 = MPI_Wtime(); // termina a contagem do tempo
  MPI_Finalize();
  if (my_rank == 0 ) {
    display( part );       // show initial step
    printf("\nTempo de execucao: %f\n\n", t2-t1); 
  }



  return 0;



}
