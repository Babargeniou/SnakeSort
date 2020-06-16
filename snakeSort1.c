//============================================================================
// Name        : snakeSort.c
// Authors     : Babar Bashir
// Version     : 1.0
// Command     : 
// Debud cmd   : mpicc snakeSort1.c -o snakeSort1
// Run		     : mpirun -n 2 snakeSort1 input.txt	-> run			
// Description : Shear sort - Parallel using MPI
//============================================================================


#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include "snake.h"

#define swap(x,y) {int t; t = x; x = y; y = t;}
#define asending(x ,y) if (x > y) swap(x , y)
#define desending(x,y) if (x < y) swap (x, y)

int cmpfunc (const void * a, const void * b);
void shearSort(int *array, int size, int my_rank, int comm_sz);
bool add_even_sort(int *array, int size);

//int recvbuff[20000];
//int sendbuff[20000];

int main(int argc, char **argv) {

    double startTime, stopTime, execution_Time;
    //int *global_Array;
    int **Matrix_A;
    int globalArraySize = 0;
    int MatrixSize = 0;
    int Row;
    int Col;
    int tag1 = 2;
    int tag2 = 1;
    int tag3 = 3;
    //int *recvbuff;
    int *sendbuff;
   // int *slaveArray;
    MPI_Status  status;
    




/************************Start MPI *******************************/

 startMPI(argc, argv);

/************************* matrix size **********************/
//readfile(inn_file,global_Array, my_rank);
 readMatrixSize(inn_file,size, my_rank);
 MatrixSize = size[0];
 RowCol = MatrixSize;
 globalArraySize = MatrixSize * MatrixSize;
 
 #ifdef DEBUG
 if (my_rank == MASTER) {
 printf("MatrixSize %d\n",size[0]);
 printf("MatrixSize %d\n",MatrixSize);
 printf("MatrixSize/RowCol %d\n",RowCol);
 printf("globalArraySize %d\n",globalArraySize);
 }
 #endif
 
/************************Reading file ****************************/
//global_Array = memory1D(500);
readfile(inn_file,global_Array, my_rank);
 
//int **Matrix_A;
Matrix_A = allocate_array(RowCol);	
createMatrix(RowCol, RowCol, global_Array, Matrix_A); 

    
#ifdef DEBUG

if (my_rank == MASTER) {
    print_matrix(Matrix_A, RowCol);

}

#endif

/* memory allocation  */
sendbuff = memory1D(globalArraySize);
//recvbuff = memory1D(globalArraySize);
//slaveArray = memory1D(globalArraySize);

/* Master Task only */

if(my_rank == MASTER){
  printf("No. of processes: %d, my rank: %d\n",comm_sz, my_rank);
  int numSlaves = comm_sz -1;
  printf("numSlaves: %d\n", numSlaves);
  
  int start = 0;
  Row = RowCol;
  Col = RowCol;
  int stop = Row*Col;
 // printf("Col %d", Col);
  int row_col_array[2]; 
  row_col_array[0]= Row;
  row_col_array[1]= Col;
  //int sendbuff[100];

  for(int dest = 1; dest < comm_sz; dest++) {
    // every processor will get row*col elements
    for(int j= start, k=0; j <stop; j++,k++){
      sendbuff[k] = *(global_Array+j);
    }

    start += (Row*Col);
		stop  += (Row*Col);
    
    // printf("message to be sent\n");
    for(int i = 0; i< globalArraySize; i++){
      printf("%d ",sendbuff[i]);
    }

   // sendbuff = memory1D(200);
    printf("\nMaster sending messages to %d\n",dest);
    MPI_Send(row_col_array, 2, MPI_INT, dest, tag1, MPI_COMM_WORLD);  /* Send rows and cols to slaves with tag 0*/
    if (MPI_Send(sendbuff, globalArraySize, MPI_INT, dest, tag2, 
                                                      MPI_COMM_WORLD) != MPI_SUCCESS) /* Send buffer to slaves with tag 1 */
    {
      printf("Error in MPI_Send when sending from rank 0.\n");
      return -1;
    }

    
  }

  /* Wait to receive results from each task */

  // wait till receives messages from all 

 // int *slaveArray;   
		//slaveArray = (int *) malloc((Row*Col)*sizeof(int));
  int *slaveArray = memory1D(16);
		
	//	int matrix[numSlaves][Row*Col];
		
//		int *merge_array;
	//	int merge_array_length = RowCol*RowCol*numSlaves;
//		merge_array = (int *) malloc(merge_array_length*sizeof(int));
	//	int start_mem_add = 0;
//		for(int i=0;i<numSlaves;i++)
	//	{
		//	printf("\nMASTER MASTER MASTER waiting MASTER MASTER MASTER %d\n ",start_mem_add);
			MPI_Recv(slaveArray, Row*Col, MPI_INT, 1, tag3, MPI_COMM_WORLD, &status);
	//		printf("\n@MASTER @MASTER @MASTER slave input received @MASTER @MASTER @MASTER\n");
		//	for(int j=0;j<(Row*Col);j++)
	//		{
	//			matrix[i][j] = slaveArray[j];
	//		}
	//		start_mem_add += (Row*Col);
//		}
	//	printf("\n");

  //  int index = 0;
  //  int start_value = 0;
  //  int end_value = Col;

   // for (int i=0; i< Row; i++) {
   //   for (int j=0; j< numSlaves; j++) {
    //    for (int k= start_value; k<end_value; k++){
      //    merge_array[index]=matrix[j][k];
      //    index = index +1;
    //    }
   //   }
    //  start_value = start_value + RowCol;
   //   end_value = end_value + RowCol;
  //  }

/* 2D merge  */    

  

} // end Master

/*  Non-Master tasks only */
else if(my_rank > MASTER){
  printf("This is slave node receive from Master\n");

  int RowColRecv[2];
  int source = MASTER;
  MPI_Recv(RowColRecv, 2, MPI_INT, source, tag1, MPI_COMM_WORLD, &status);
   Row = RowColRecv[0];
   Col = RowColRecv[1];
   //printf("globalArraySize: %d\n", Row*Col);
   int *recvbuff = memory1D(Row*Col);
  //int recvbuff[100];
  //int globArraySize = Row * Col;
 // printf("globalArraySize: %d\n", globArraySize);
  MPI_Recv(recvbuff, Row*Col, MPI_INT, source, tag2, MPI_COMM_WORLD, &status);
  
  //MPI_Barrier( MPI_COMM_WORLD );
  int recv_count;
    printf("receive message\n");
    for(int i = 0; i< Row*Col; i++){
      printf("check1: %d ",recvbuff[i]);
    }
    printf("\n");
  MPI_Get_count(&status, MPI_INT, &recv_count);

  printf("recv_count: %d\n", recv_count);

  /* shear sort */
 // shearSort(recvBuff, RowCol,  my_rank,  comm_sz);
  int desti = MASTER;
  printf("\nSENDING MESSAGE BACK TO MASTER\n");	//TAG 2
  if (MPI_Send(recvbuff, recv_count, MPI_INT, desti, tag3, MPI_COMM_WORLD) != MPI_SUCCESS)
  {
    printf("Sending back error");
  }
  
  startTime = MPI_Wtime();

  printf("receive message\n");
    for(int i = 0; i< Row*Col; i++){
      printf("check2: %d ",recvbuff[i]);
    }
    printf("\n");

  stopTime = MPI_Wtime();

  execution_Time = stopTime - startTime;

  printf("Execution time: %f\n", execution_Time);
  

} // slaves rank





free_array(Matrix_A, RowCol);
//free(global_Array);
//free(recvbuff);
free(sendbuff);
//free(slaveArray);

done();

 
return 0;
}  // end main



int cmpfunc (const void * a, const void * b) {
   return ( *(int*)b - *(int*)a );
}


void shearSort(int *array, int size, int my_rank, int comm_sz) {

	// memory for local array

	//int d = ceil(log10(size)/log10(2));
	int d = ceil(log2(size));
  int rows = size;
  int cols = size;
  bool sorted = false;
  bool endsort = false;
  

	for(int l = 1; l < d ; l++ ) {
    
    for (int i = 0; i < rows; i+=2){  //  even sort
       
      int *evenRowArray; 
      evenRowArray = memory1D(size);

      for (int j = 0; j < cols; j++ ){
      /*build an array which is need to be odd-even sorted.*/
      evenRowArray[j] = *(array + (i*cols)+j);
      }

      // odd_even sort function.

      for (int j = 0; j<cols; j++){
        /* put values back in orginal array */
        *(array +(i*cols)+j) = evenRowArray[j];
      }


		}
    
    for(int i = 1; i < rows; i= i+2){ // odd sort
     
      int *oddRowArray; 
      oddRowArray = memory1D(size);

      for (int j = 0; j < cols; j++ ){
      /*build an array which is need to be odd-even sorted.*/
      oddRowArray[j] = *(array + (i*cols)+j);
      }

      // odd_even sort function.
      
      //Reverse array
			for(int i=0,j=cols-1;i<cols/2;i++,j--){
				int tmp = oddRowArray[i];
				oddRowArray[i] = oddRowArray[j];
				oddRowArray[j] = tmp;
			}

      for (int j = 0; j< cols; j++){
        /* put values back in orginal array */
        *(array +(i*cols)+j) = oddRowArray[j];
      }
			

  	}

    for (int i = 0; i< cols; i++) {
      // column sorting
      int * columnArray;
      columnArray = memory1D(size);

      for(int j = 0; j< rows; j++) {

        columnArray[j] = *(array+i+(j*cols));

      }

       // odd_even sort function.

       for(int j = 0; j< rows; j++) {
         // put values back in orginal column
         *(array+i+(j*cols)) = columnArray[j];
       }

    }
  }
}


bool add_even_sort(int *array, int size) 
{

  for (int i = 0; i<size; i++){
    // print array here
  }

  bool sorted = false;
  bool endsort = true;  

  while(!sorted) {
    sorted = true;
    for(int i = 1; i< size-1; i+=2) {
      if(*(array+i) > *(array+i+1))
      {   
      // odd sorting
      int temp = *(array+i);
      *(array+i) = *(array+i +1);
      *(array+i+1) = temp;
      sorted = false;
      endsort = false;
      }
    }

    for (int i = 0; i < size-1; i+=2){
      if(*(array+i) > *(array+i+1))
      {
      // even sorting
        int temp = *(array+i);
        *(array+i) = *(array+i +1);
        *(array+i+1) = temp;
        sorted = false;
        endsort = false;
      }
    }
  }
  return endsort;
}