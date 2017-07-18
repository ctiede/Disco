#include "paul.h"

double get_dV( double *, double * );
void printTracerCoords( struct domain * )i;
void addTracers( struct tracerList * );
int getN0( int , int , int );

void setTracerParams( struct domain * theDomain){

   int num_tracers = theDomain->theParList.num_tracers; 
   int Nmc = num_tracers;
   theDomain->Nmc = Nmc;
   
   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;
   int *Np = theDomain->Np;
   int jk;
   for( jk=0; jk<Nr*Nz; ++jk )
      tr_tot += Np[jk]*Nmc;
   theDomain->Ntr = tr_tot;
}

void initMCtracers( struct domain *theDomain ){

   int Nmc = theDomain->Nmc;
   int Ntr = theDomain->Ntr;
   int Nr = theDomain->Nr;
   int Nz = theDoamin->Nz;
   int *Np = theDomain->Np;
   int rank = theDomain->rank;
   int size = theDomain->size;

   int l;
   int id = 0;
   int count = 0;
   for( l=0; l<size; ++l ){
      MPI_Barrier();
      if( l==rank ){
         id = count;
         count = Ntr;
         MPI_Allreduce(MPI_IN_PLACE, &count, 1, MPI_INT, MPI_SUM, theDomain->theComm);
      }
   }
   
   //function that takes ID and theCells and gives the cells their tracers ->or just do it here?
   int n,i,j,k;
   for( j=0; j<Nr; ++j ){
      for( k=0; k<Nz; ++k ){
         jk = j + Nr*k;
         for( i=0; i<Np[jk]; ++i ){
            struct cell *c = &(theCells[jk][i]);
            for( n=0; n<Nmc; ++n )
               addTracer( c, id );
         }
      }
   }
}

int check_phi(double phi, double phip, double dphi, double phi_max){

   double phic = phi - phip;
   while( phic >  phi_max/2.0 ){ phic -= phi_max; }
   while( phic < -phi_max/2.0 ){ phic += phi_max; }
   if( -dphi<phic && phic<0 ){ return 1; }
   else{ return 0; }

}

int check_in_cell(struct tracer *tr, double *xp, double *xm, double phi_max){

   double r   = tr->R;
   double phi = tr->Phi;
   double z   = tr->Z;

   if( xm[0] < r && r < xp[0] ){
	if( xm[2] < z && z < xp[2]){
	     if( check_phi(phi, xp[1], xp[1]-xm[1], phi_max) ){
   		return 1;
	     }
	}
   }
   return 0;
}

void test_cell_vel( struct tracer *tr, struct cell *c ){

   int check1=0, check2=0, check3=0;
   if( isnan(c->prim[URR]) ){
	   tr->Vr 	  = 0;
	   tr->Omega = 0;
   	tr->Vz    = 0;
	   tr->Type  = 3;
   	check1 = 1;
   }
   if( isnan(c->prim[UPP]) ){
	   tr->Vr    = 0;
      tr->Omega = 0;
      tr->Vz    = 0;
      tr->Type  = 4;
   	check2 = 1;
   }
   if( isnan(c->prim[UZZ]) ){
	   tr->Vr    = 0;
      tr->Omega = 0;
      tr->Vz    = 0;
      tr->Type  = 5;
	   check3 = 1;
   }

  if( check1==1 && check2==1 && check3==1){ tr->Type = 6; }

}


void get_local_vel(struct tracer *tr, struct cell *c){

   double vr, om, vz;
   int type = 1;

  if( c != NULL ){
      vr = c->prim[URR];
      om = c->prim[UPP];
      vz = c->prim[UZZ];
   } else{
      vr = 0.0;
      om = 0.0;
      vz = 0.0;
      type = 2;
   }
   tr->Vr    = vr;
   tr->Omega = om;
   tr->Vz    = vz;
   //tr->Type  = type;
   if( c != NULL ){
     //test_cell_vel( tr, c );
   }
}

struct cell * get_tracer_cell(struct domain *theDomain, struct tracer *tr){

   struct cell **theCells = theDomain->theCells;
   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;
   int *Np = theDomain->Np;
   double phi_max = theDomain->phi_max;
   double *r_jph = theDomain->r_jph;
   double *z_kph = theDomain->z_kph;

   int i,j,k;
   for( j=0; j<Nr; ++j){
	for( k=0; k<Nz; ++k){
	   int jk = j+Nr*k;
	   double rm = r_jph[j-1];
	   double rp = r_jph[j];
	   double zm = z_kph[k-1];
	   double zp = z_kph[k];
	   for( i=0; i<Np[jk]; ++i){
	   	struct cell *c = &(theCells[jk][i]);
	   	double phip = c->piph;
		   double phim = phip - c->dphi;
		   double xp[3] = {rp, phip, zp};
		   double xm[3] = {rm, phim, zm};
		   if( check_in_cell(tr, xp, xm, phi_max) ){
			   return c;
		}
	   }
	}
   }
   return NULL;

}


void moveTracers(struct domain *theDomain, struct tracer *tr, double dt){

   double rmin = theDomain->theParList.rmin;
   double rmax = theDomain->theParList.rmax;
   double zmin = theDomain->theParList.zmin;
   double zmax = theDomain->theParList.zmax;

   double r = tr->R;
   double phi = tr->Phi;
   double z = tr->Z;

   r += tr->Vr*dt;
   if( r > rmax ) r = 0.0/0.0;
   if( r < rmin ) r = 0.0/0.0;

   z += tr->Vz*dt;
   if( z > zmax ) z = 0.0/0.0;
   if( z < zmin ) z = 0.0/0.0;

   tr->R = r;
   tr->Z = z;

   phi += tr->Omega*dt;
   tr->Phi = phi;
}


void tracer_RK_copy( struct tracer * tr ){
   tr->RK_r     = tr->R;
   tr->RK_phi   = tr->Phi;
   tr->RK_z	= tr->Z;
   tr->RK_vr	= tr->Vr;
   tr->RK_omega = tr->Omega;
   tr->RK_vz	= tr->Vz;
}

void tracer_RK_adjust( struct tracer * tr , double RK ){
   tr->R     = (1.-RK)*tr->R     + RK*tr->RK_r;
   tr->Phi   = (1.-RK)*tr->Phi   + RK*tr->RK_phi;
   tr->Z     = (1.-RK)*tr->Z     + RK*tr->RK_z;
   tr->Vr    = (1.-RK)*tr->Vr    + RK*tr->RK_vr;
   tr->Omega = (1.-RK)*tr->Omega + RK*tr->RK_omega;
   tr->Vz    = (1.-RK)*tr->Vz    + RK*tr->RK_vz;
}


void updateTracers(struct domain *theDomain, double dt){

   //printf("Trying to update tracers!");
   struct tracer *tr = theDomain->theTracers->head;
   while( tr!=NULL ){
   	struct cell   *c  = get_tracer_cell( theDomain , tr );
	   get_local_vel( tr , c );
      moveTracers( theDomain , tr , dt );
      tr = tr->next;
   }

}
