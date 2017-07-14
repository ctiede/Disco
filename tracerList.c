#include "paul.h"


void addTracer( struct tracerList *theList ){
   //adds tracer to the front of the list
   struct tracer *tr = (struct tracer *) malloc( sizeof(struct tracer) );

   tr->next = theList->head;

   theList->head = tr;
   theList->size += 1;
}

int getListSize( struct tracerList *theList ){

   int size = 0;
   struct tracer *tr = theList->head;
   while( tr!=NULL ){
      size++;
      tr = tr->next;
   }
  return size;
}

void rmTracers( struct tracerList *theList ){

  struct tracer *tr = theList->head;
  struct tracer *prev = NULL;
  while( tr!=NULL ){
      if( tr->rmFlag ){
         if( prev ){ //if prev has a valid value
            prev->next = tr->next;
         }
         else{ //if prev is NULL, need to repoint head
            theList->head = tr->next;
         }
         free(tr);
      }
      prev = tr;
      tr = tr->next;
   }
}

void init_tracerList( struct domain *theDomain ){

   int Ntr = theDomain->Ntr;
   struct tracerList *theList = theDomain->theTracers;
   theList->head = NULL;
   theList->size = 0;

   int n;
   printf("Trying to create list of %d tracers\n", Ntr);
   for( n=0; n<Ntr; ++n){     //create linked list of Ntr tracers
      addTracer( theList );
   }

   //printf("Trying to get list size (%d were added)\n", theList->size);
   //int size = getListSize( theList );
   //printf("List Size: %d\n", size);

}

void printTracerCoords( struct domain *theDomain ){

   struct tracer *tr = theDomain->theTracers->head;
   int type;
   double r, z, phi;
   printf("Type, (r, phi, z)\n");
   while ( tr!=NULL ){
      type = tr->Type;
      r = tr->R; z = tr->Z; phi = tr->Phi;
      printf("%d,\t (%4.4f,\t %4.4f,\t %4.4f)\n", type, r, phi, z);
      tr = tr->next;
   }
}
