#include "paul.h"

int getN0( int, int, int );

void addTracer( *theList ){

   struct tracer *temp = theList->head;
   struct tracer tr = {0};
   theList->head = tr;
   tr.next = temp;
}

void init_tracerList( struct domain *theDomain ){

   int Ntr = theDomain->Ntr;
   struct tracerList *theList = theDomain->theTracers;
   theList->head = NULL;

   int n;
   for( n=0; n<Ntr; ++n){     //create linked list of Ntr tracers
      addTracer( theList );
   }

}

int getListSize( *theList ){

   int size = 0;
   struct tracer *tr = theList->head;
   while( tmp != NULL){
      size++;
      tr = tr->next;
   }
  return size;
}

