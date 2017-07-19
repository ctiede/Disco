#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "paul.h"

struct tr_lite{

   int type;
   double x[3];
   double v[3];
};

void generate_mpi_tr( MPI_Datatype *tr_mpi ){

   struct tracer_lite test;
   int count = 3;
   int blocksize[]      = { 1, 3, 3 };
   MPI_Datatype types[] = { MPI_Int, MPI_DOUBLE, MPI_DOUBLE };
   MPI_Aint offsets[3];
   
   offsets[0] = (char *)&(test.type) - (char *)(&test);
   offsets[1] = (char *)&(test.x)    - (char *)(&test);
   offsets[2] = (char *)&(test.v)    - (char *)(&test);

   MPI_Type_create_struct( count, blocksize, offsets, types, tr_mpi );
   MPI_Type_commit( tr_mpi );
}

void copy_tr_to_lite( struct tracer *tr, struct tr_lite *trl ){

   memcpy( trl->x, tr->x, 3*sizeof(double) );
   memcpy( trl->v, tr->v, 3*sizeof(double) );
   trl->type = tr->type;

}

void copy_lite_to_tr( struct tr_lite *trl, struct tracer *tr ){

   memcpy( tr->x, trl->x, 3*sizeof(double) );
   memcpy( tr->v, trl->v, 3*sizeof(double) );
   tr->type = trl->type;

}

void generate_sendbuffer_tr( struct domain *theDomain, int numL, int numR,  struct tr_lite *trL, struct tr_lite *trR, dim ){
   
   int iL = 0;
   int iR = 0;
   struct tracerList *theList = theDomain->theTracers;
   while( tr!=NULL ){
      if( tr->rmFlag==1 ){
         copy_tr_to_lite( tr, trL+iL );
         iL++
      }
      if( tr->rmFlag==2 ){
         copy_tr_to_lite( tr, trR+iR );
         iR++
      }
   }

}

void generate_tr_number( struct domain *theDomain, int *numL, int *numR, int dim ){

   double r0   = theDomain->r0;
   double z0   = theDomain->z0;
   double delr = theDomain->delr;
   double delz = theDomain->delz;
   double x0   = { r0  , z0  };
   double delx = { delr, delz};

   int iL = 0;
   int iR = 0;
   double coord;
   struct tracer *tr = theDomain->theTracers->head;
   while( tr!=NULL ){
      if( dim==0 ) 
         coord = tr->R;
      if( dim==1 ) 
         coord = tr->Z;
      if( coord < x[dim] ){
         tr->rmFlag = 1;
         iL++;         
      }
      if( x[dim]+delx[dim] < coord ){
         tr->rmFlag = 2;
         iR++;
      }
      tr = tr->next;
   }
   *numL = iL;
   *numR = iR;

}

void addTracer( struct tracerList * );

void exchangeTracers( struct *theDomain, int dim ){
   
   MPI_Datatype tr_mpi = {0};
   generate_mpi_tr( &tr_mpi );

   MPI_Comm grid_comm = theDomain->theComm;
   int *left_rank = theDomain->left_rank;
   int *right_rank = theDomain->right_rank;

   int n;
   int tag = 0;
   MPI_Status status;
   struct tracerList theList = theDomain->theTracers;

   int numL, numR;
   int send_sizeL = 0;
   int send_sizeR = 0;
   int recv_sizeL = 0;
   int recv_sizeR = 0;

//------Count Tracers to Send----------------------------
   generate_tr_number( theDomain, &numL, &numR, dim );
   send_sizeL = numL;
   send_sizeR = numR;

//------Tell Neighbors how many to expect----------------
   MPI_Sendrecv( &send_sizeL, 1, MPI_INT,  left_rank[dim], tag,
                 &recv_sizeR, 1, MPI_INT, right_rank[dim], tag,
                 grid_comm, &status );
   MPI_Sendrecv( &send_sizeR, 1, MPI_INT, right_rank[dim], tag+1,
                 &recv_sizeL, 1, MPI_INT,  left_rank[dim], tag+1,
                 grid_comm, &status );

//------Build list of tracers to send--------------------
   struct tr_lite trL_send[send_sizeL];
   struct tr_lite trR_send[send_sizeR];
   struct tr_lite trL_recv[recv_sizeL];
   struct tr_lite trR_recv[recv_sizeR];

//------Build list of tracers to send---------------------
   generate_sendbuffer_tr( theDomain, numL, numR, trL_send, trR_send, dim );

//------Send lists of tracers----------------------------- 
   MPI_Sendrecv( trL_send, send_sizeL, tr_mpi,  left_rank[dim], tag+2,
                 trR_recv, recv_sizeR, tr_mpi, right_rank[dim], tag+2,
                 grid_comm, &status );
   MPI_Sendrecv( trR_send, send_sizeR, tr_mpi, right_rank[dim], tag+3,
                 trL_recv, recv_sizeL, tr_mpi,  left_rank[dim], tag+3,
                 grid_comm, &status );
//------Add/Clear tracers from theDomain
   for( n=0; n<recv_sizeR; ++n){
      addTracer( theDomain->theTracers );
      //function to add tr_lite info to the tracer that has just been added
   }

   MPI_Type_free( &tr_mpi );
}
