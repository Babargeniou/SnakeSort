//============================================================================
// Name        : snakeSort.c
// Authors     : Babar Bashir
// Version     : 1.0
// Command     : 
// Debud cmd   : g++ -g -fopenmp q4.cpp -o out
// Run		   : mpirun -n 2 snakeSort input.txt	-> run			
// Description : Shear sort - Parallel using MPI
//============================================================================


#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "snake.h"

#define swap(x,y) {int t; t = x; x = y; y = t;}
#define asending(x ,y) if (x > y) swap(x , y)
#define desending(x,y) if (x < y) swap (x, y)

int cmpfunc (const void * a, const void * b);

int main(int argc, char **argv) {

int sum = 0;
int rec_buf[100];          // buffer where the received data should be stored

/************************Start MPI *******************************/

 startMPI(argc, argv);

/************************* matrix size **********************/
// readMatrixSize(inn_file,size);
 //printf("size %d\n",size[0]);
 
// int SIZE = size[0];
 //printf("size %d\n",SIZE);
 
 
 //int *Array;
 //Array = (int*)malloc((SIZE*SIZE) * sizeof(int));
  // if (Array == NULL) {
  //     printf("Memory not allocated for Matrixx \n");
  // } 
 //  else
 //  {
 //      printf("Memory allocated successfully for Matrixx\n");
 //  }

/************************Reading file ****************************/
 readfile(inn_file,Array);
 //RowCol =(int) x[-1];
 readMatrixSize(inn_file,size);
 RowCol = size[0];

int **Matrix_A;
Matrix_A = allocate_array(RowCol);	


createMatrix(RowCol, RowCol, Array, Matrix_A); 


chunk  = RowCol/comm_sz;                 /* Number of intervals per processor   */
//if (my_rank == MASTER) {
//printf("localSize: %d\n",chunk);
//}
if (my_rank == MASTER) {
print_matrix(Matrix_A, RowCol);
}


int rem = (RowCol*RowCol)%comm_sz; // elements remaining after division among processes

    int *sendcounts = malloc(sizeof(int)*comm_sz);
    int *displs = malloc(sizeof(int)*comm_sz);

    // calculate send counts and displacements
    for (int i = 0; i < comm_sz; i++) {
        sendcounts[i] = (RowCol*RowCol)/comm_sz;
        if (rem > 0) {
            sendcounts[i]++;
            rem--;
        }

        displs[i] = sum;
        sum += sendcounts[i];
    }

    // print calculated send counts and displacements for each process
    if (my_rank == MASTER) {
        for (int i = 0; i < comm_sz; i++) {
            printf("sendcounts[%d] = %d\tdispls[%d] = %d\n", i, sendcounts[i], i, displs[i]);
        }
    }

    // divide the data among processes as described by sendcounts and displs
    MPI_Scatterv(&Array, sendcounts, displs, MPI_INT, &rec_buf, 100, MPI_INT, 0, MPI_COMM_WORLD);

    // print what each process received
    printf("%d: ", my_rank);
    for (int i = 0; i < sendcounts[my_rank]; i++) {
        //printf("%d\t", rec_buf[i]);
        qsort(rec_buf, 8, sizeof(int), cmpfunc);
       // printf("\n");
        printf("%d\t", rec_buf[i]);

    }
    printf("\n");

    int **Matrix_A_Sort;
Matrix_A_Sort = allocate_array(RowCol);	
    createMatrix(RowCol, RowCol, rec_buf, Matrix_A_Sort); 

 //   if (my_rank == 1) {
//7print_matrix(Matrix_A_Sort, RowCol);
//}

if (my_rank == 0) {
print_matrix(Matrix_A_Sort, RowCol);
}



/***** calculate log base 2 ****************/
/****** Number of steps d +1 **************/
int d = ceil(log10(RowCol)/log10(2));

//if (my_rank == MASTER) {
  //  printf("d: %d\n", d);
//}


 /* free array   */
 free_array(Matrix_A, RowCol);
 //free(Array);


 done();
 return 0;
}  // end main



int cmpfunc (const void * a, const void * b) {
   return ( *(int*)b - *(int*)a );
}