#include "paul.h"

int getN0( int, int, int );

void init_tracerList( struct domain *theDomain ){

//--------Put this into a function -> give theDomain vars N0r, N0z, delR, delZ------------------
   int Num_R = theDomain->theParList.Num_R;
   int Num_Z = theDomain->theParList.Num_Z;
   int phi_max = theDomain->theParList.phimax;
   int *dim_rank = theDomain->dim_rank;
   int *dim_size = theDomain->dim_size;

   int N0r = getN0( dim_rank[0], dim_size[0], Num_R );
   int N1r = getN0( dim_rank[0]+1, dim_size[0], Num_R);
   int Nr = N1r - N0r;

   int N0z = getN0( dim_rank[1], dim_size[1], Num_Z );
   int N1z = getN0( dim_rank[1]+1, dim_size[1], Num_Z);
   int Nz = N1z - N0zi;

   double dr = 1.0/(double)Num_R;
   double r0 = (double)N0r*dr;
   double delR = Nr*dr;

   double dz = 1./(double)Num_Z;
   double z0 = (double)N0z*dz;
   double delZ = Nz*dz;
//-----------Save these things in gridsetup.c------------------------------------------------

   srand(theDomain->rank);
   rand();

   int Ntr = theDomain->Ntr;
   struct tracerList *theList = theDomain->theTracers;

   //struct tracer tr[Ntr];           //array of tracers for initializing
   //struct tracer *temp = &tr[0];
   theList->head = NULL;

   int n;
   for( n=0; n<Ntr; ++n){     //create linked list of Ntr tracers
      addTracer( theList );
   }
   
}

void addTracer( *theList ){

   struct tracer *temp = theList->head;
   struct tracer tr = {0};
   theList->head = tr;
   tr.next = temp;
}
