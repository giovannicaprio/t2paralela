/*  
 Game of life
 D. Thiebaut
 This is the MPI version of GameOfLife.c, for MPI, for 2 Processes,
 or tasks.
 This version works only for 2 tasks.

 It hasn't been optimized.  Each task works on one half of the dish
 arrays where the generations evolve, but actually maintain a full
 array in memory.
 
 The two tasks communicate with each other after each generation 
 and exchange 2 rows, each.  

 This code works in console mode, displaying successive generations
 of the game of life on the screen, and clearing the screen between
 each one.

 The initial pattern is defined in the array PATTERN.

 To compile and run, several different variants:

 mpicc.mpich -o GameOfLifeMPI2 GameOfLifeMPI2.c
 mpiexec.mpich -n 2 ./GameOfLifeMPI2
 
 or 

 mpicc -o GameOfLifeMPI2 GameOfLifeMPI2.c
 mpirun -np 2 ./GameOfLifeMPI2

*/
#include <stdio.h>
#include <strings.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include "mpi.h"

#define NUMBERROWS 28
#define esc 27
#define cls() printf("%c[2J",esc)
#define pos(row,col) printf("%c[%d;%dH",esc,row,col)


char  *DISH0[ NUMBERROWS ];
char  *DISH1[ NUMBERROWS ];
char  *PATTERN[NUMBERROWS] = {
  "                                                                                  ",
  "   #                                                                              ",
  " # #                                            ###                               ",
  "  ##                                                                              ",
  "                                                                                  ",
  "                                                      #                           ",
  "                                                    # #                           ",
  "                                                     ##                           ",
  "                                                                                  ",
  "                                                                                  ",
  "                                                                                  ",
  "                                                                                  ",
  "             #                                                                    ",
  "           # #                                                                    ",
  "            ##                                                                    ",
  "                                                                                  ",
  "                                                                                  ",
  "                                                                                  ",
  "                                                                                  ",
  "                                                                                  ",
  "                                                                                  ",
  "                                                                                  ",
  "                                                                                  ",
  "                                                                                  ",
  "                                                                                  ",
  "                                                                                  ",
  "                                                                                  ",
  "                                                                                  " 
};
int ROWSIZE = strlen( "                                                                                  ") + 1;

//------------------------------- prototypes --------------------------------
void life( char**, char**, int );
void initDishes( int  );
void print( char **, int );

// --------------------------------------------------------------------------
// initDishes
// inits the dishes (current and future)
void initDishes( int rank ) {
  int i;

  //--- initialize other dish with spaces.  Make it same dimension as DISH0. ---
  for (i = 0; i< NUMBERROWS; i++ )  {
    DISH0[i] = (char *) malloc( ( strlen( PATTERN[0] ) + 1 ) * sizeof( char ) );
    strcpy( DISH0[i], PATTERN[i] );

    DISH1[i] = (char *) malloc( (strlen( DISH0[0] )+1) * sizeof( char )  );
    strcpy( DISH1[i], PATTERN[i] );
  }
}


// --------------------------------------------------------------------------
// initDishes2
// inits the dishes (current and future)
// (Buggy: attempts to declare only 1 half of the array, plus boundary
//  rows, depending on rank.  Needs a bit more debugging...)
void initDishes2( int rank ) {
  int i;

  // init to null all entries.  This way we'll be
  // able to tell if a row belongs to us or not.
  for ( i=0; i<NUMBERROWS; i++ ) {
    DISH0[i] = NULL;
    DISH1[i] = NULL;
  }

  //--- Init RANK 0 rows ---
  if ( rank == 0 ) {
    //--- initialize dishes with lower half of pattern ---
    for (i = 0; i< NUMBERROWS/2+1; i++ )  {
    
      DISH0[i] = (char *) malloc( (strlen( PATTERN[0] ) + 1 ) * sizeof( char ) );
      strcpy( DISH0[i], PATTERN[i] );

      DISH1[i] = (char *) malloc( (strlen( DISH0[0] )+1) * sizeof( char )  );
      strcpy( DISH1[i], DISH0[i] );
    }

    //--- initialize top row of dishes, as they are neighbors of Row 0 ---
    DISH0[NUMBERROWS-1] = (char *) malloc( (strlen( PATTERN[0] ) + 1 ) * sizeof( char ) );
    strcpy( DISH0[NUMBERROWS-1], PATTERN[NUMBERROWS-1] );
    DISH1[NUMBERROWS-1] = (char *) malloc( (strlen( PATTERN[0] ) + 1 ) * sizeof( char ) );
    strcpy( DISH1[NUMBERROWS-1], PATTERN[NUMBERROWS-1] );
  }

  //--- Init RANK 1 rows ---
  if ( rank == 1 ) {
    //--- initialize dishes with upper half of pattern ---                                                                               
    for (i = NUMBERROWS/2-1; i< NUMBERROWS; i++ )  {

      DISH0[i] = (char *) malloc( (strlen( PATTERN[0] ) + 1 ) * sizeof( char ) );
      strcpy( DISH0[i], PATTERN[i] );

      DISH1[i] = (char *) malloc( (strlen( DISH0[0] )+1) * sizeof( char )  );
      strcpy( DISH1[i], DISH0[i] );
    }

    //--- initialize bottom row of dishes, as they are neighbors of top row  ---                                                       
    DISH0[0] = (char *) malloc( (strlen( PATTERN[0] ) + 1 ) * sizeof( char ) );
    strcpy( DISH0[0], PATTERN[0] );
    DISH1[0] = (char *) malloc( (strlen( PATTERN[0] ) + 1 ) * sizeof( char ) );
    strcpy( DISH1[0], PATTERN[01] );
  }
}



// --------------------------------------------------------------------------
// print
void print( char* dish[], int rank ) {
  int i;

  if ( rank == 0 ) {
    //--- display lower half only ---
    for (i=0; i<NUMBERROWS/2; i++ ) {
      if ( dish[i] == NULL ) continue;
      pos( i, 0 );
      printf( "%s\n", dish[i] );
    }
  }

  if ( rank == 1 ) {
    //--- display upper half only ---
    for (i=NUMBERROWS/2; i<NUMBERROWS; i++ ) {
      if ( dish[i] == NULL ) continue;
      pos( i, 0 );
      printf( "%s\n", dish[i] );
    }

  }
}

// --------------------------------------------------------------------------
void check( char** dish, char** future ) {
  int i, j, k, l;
  l = sizeof( dish )/sizeof( dish[0] );
  printf( "length of dish = %d\n", l );

  for ( i=0; i<l; i++ ) {
    k = strlen( dish[i] );
    printf( "%d %s\n", k, dish[i] );
  }
  printf( "\n\n" );

  l = sizeof( future )/sizeof( future[0] );
  printf( "length of future = %d\n", l );

  for ( i=0; i<l; i++ ) {
    k = strlen( future[i] );
    printf( "%d %s\n", k, future[i] );
  }
  printf( "\n\n" );

}

// --------------------------------------------------------------------------
void  life( char** dish, char** newGen, int rank ) {
  /*
   * Given an array of string representing the current population of cells
   * in a petri dish, computes the new generation of cells according to
   * the rules of the game. A new array of strings is returned.
   */
  int i, j, row;
  int rowLength = strlen( dish[0] );
  int dishLength = NUMBERROWS;
  
  int lowerRow, upperRow;

  //--- slice the array into two halves.  Rank 0 is lower half, ---
  //--  Rank 1 is upper half.                                   ---
  if ( rank == 0 ) {
    lowerRow = 0;
    upperRow = NUMBERROWS/2;
  }
  if ( rank == 1 ) {
    lowerRow = NUMBERROWS/2;
    upperRow = NUMBERROWS;
  }
  
  for (row = lowerRow; row < upperRow; row++) {// each row
    
    if ( dish[row] == NULL )
      continue;

    for ( i = 0; i < rowLength; i++) { // each char in the
                                       // row

      int r, j, neighbors = 0;
      char current = dish[row][i];

      // loop in a block that is 3x3 around the current cell
      // and count the number of '#' cells.
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
          if (dish[realr][realj] == '#')
            neighbors++;
        }
      }

      if (current == '#') {
        if (neighbors < 2 || neighbors > 3)
          newGen[row][i] =  ' ';
        else
          newGen[row][i] = '#';
      }

      if (current == ' ') {
        if (neighbors == 3)
          newGen[row][i] = '#';
        else
          newGen[row][i] = ' ';
      }
    }
  }
}


// --------------------------------------------------------------------------
int main( int argc, char* argv[] ) {
  int gens = 3000;      // # of generations
  int i;
  char **dish, **future, **temp;

  //--- MPI Variables ---
  int noTasks = 0;
  int rank = 0;
  MPI_Status status;   // required variable for receive routines

  //--- initialize MPI ---
  MPI_Init( &argc, &argv );

  //--- get number of tasks, and make sure it's 2 ---
  MPI_Comm_size( MPI_COMM_WORLD, &noTasks );
  if ( noTasks != 2 ) {
    printf( "Number of Processes/Tasks must be 2.  Number = %d\n\n", noTasks );
    MPI_Finalize();
    return 1;
  }

  //--- get rank ---
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );

  //--- init the dishes as half of the original problem ---
  initDishes( rank );

  dish   = DISH0;
  future = DISH1;

  //check( dish, future );

  //--- clear screen ---
  cls();

  print( dish, rank );          // # first generation, in petri dish


  // iterate over all generations
  for ( i = 0; i < gens; i++) {
    
    pos( 33+rank, 0 );
    printf( "Rank %d: Generation %d\n", rank, i );

    // apply the rules of life to the current population and 
    // generate the next generation.
    life( dish, future, rank );

    // display the new generation
    //print( dish, rank );

    // add a bit of a delay to better see the visualization
    // remove this part to get full timing.
    //if (rank == 0 ) sleep( 1 );

    if (rank==0 ) {
      //          buffer                #items   item-size src/dest tag   world  
      MPI_Send( future[    0         ], ROWSIZE, MPI_CHAR,    1,     0,  MPI_COMM_WORLD );
      MPI_Send( future[NUMBERROWS/2-1], ROWSIZE, MPI_CHAR,    1,     0,  MPI_COMM_WORLD );
      MPI_Recv( future[NUMBERROWS-1],   ROWSIZE, MPI_CHAR,    1,     0,  MPI_COMM_WORLD, &status );
      MPI_Recv( future[NUMBERROWS/2],   ROWSIZE, MPI_CHAR,    1,     0,  MPI_COMM_WORLD, &status );
    }
    if (rank==1 ) {
      MPI_Recv( future[    0         ], ROWSIZE, MPI_CHAR,    0,     0,  MPI_COMM_WORLD, &status );
      MPI_Recv( future[NUMBERROWS/2-1], ROWSIZE, MPI_CHAR,    0,     0,  MPI_COMM_WORLD, &status );
      MPI_Send( future[NUMBERROWS-1],   ROWSIZE, MPI_CHAR,    0,     0,  MPI_COMM_WORLD );
      MPI_Send( future[NUMBERROWS/2],   ROWSIZE, MPI_CHAR,    0,     0,  MPI_COMM_WORLD );
    }
    // copy future to dish
    temp = dish;
    dish = future;
    future = temp;
  }

  //--- display the last generation ---
  print(dish, rank);

  //--- close MPI ---
  pos( 30+rank, 0 );
  printf( "Process %d done.  Exiting\n\n", rank );
  MPI_Finalize();

  return 0;
}