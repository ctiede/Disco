#include "paul.h"

void setTracerParams( struct domain * theDomain){

     int num_tracers = theDomain->theParList.num_tracers;
     printf("Total Tracers: %d\n", num_tracers);
     int size = theDomain->size;
     printf("Size: %d\n", size);
     theDomain->Ntr = num_tracers/size;   //crude round-down for now
     printf("Local Tracers: %d\n", num_tracers/size);
     //theDomain->Ntr = 2000;
}

double getRandIn( double xmin, double dx){

  return xmin + ( (double)rand()/(double)RAND_MAX )*dx;
}


void printTracerCoords( struct domain * );

void initTracers_Rand( struct domain *theDomain ){  //randomly init tracers in serial

   struct param_list theParamList = theDomain->theParList;
   double rmin = theParamList.rmin;
   double rmax = theParamList.rmax;
   double dr   = rmax - rmin;
   double zmin = theParamList.zmin;
   double zmax = theParamList.zmax;
   double dz   = zmax - zmin;
   double phimax = theParamList.phimax;

   printf("Domain:\n r0: %f\n delr: %f\n z0: %f\n  delz: %f\n", rmin,dr,zmin,dz);

   srand(theDomain->rank);
   rand();

   struct tracerList *theList = theDomain->theTracers;
   struct tracer *tr = theList->head;
   while( tr!=NULL ){
     double r = getRandIn( rmin, dr );
     double z = getRandIn( zmin, dz );
     double phi = getRandIn( 0.0, phimax );
     //printf("Generated (r,phi,z): %f, %f, %f\n", r, phi, z);
     tr->R = r; tr->Z = z; tr->Phi = phi;
     tr->Type = 0; tr->rmFlag = 0;
     tr = tr->next;
   }
   //printTracerCoords( theDomain );
}


int getN0( int , int , int );

void initializeTracers( struct domain *theDomain ){

//----Put this into a function? Give domain vars N0r, N0z, delR, delZ?------
   int Num_R = theDomain->theParList.Num_R;
   int Num_Z = theDomain->theParList.Num_Z;
   int phi_max = theDomain->theParList.phimax;
   int *dim_rank = theDomain->dim_rank;
   int *dim_size = theDomain->dim_size;

   double rmin = theDomain->theParList.rmin;
   double rmax = theDomain->theParList.rmax;
   double zmin = theDomain->theParList.zmin;
   double zmax = theDomain->theParList.zmax;

   int N0r = getN0( dim_rank[0], dim_size[0], Num_R );
   int N1r = getN0( dim_rank[0]+1, dim_size[0], Num_R);
   int Nr = N1r - N0r;

   int N0z = getN0( dim_rank[1], dim_size[1], Num_Z );
   int N1z = getN0( dim_rank[1]+1, dim_size[1], Num_Z);
   int Nz = N1z - N0z;

   double dr = (rmax-rmin)/(double)Num_R;
   double r0 = rmin + (double)N0r*dr;
   double delr = Nr*dr;

   double dz = (zmax-zmin)/(double)Num_Z;
   double z0 = zmin + (double)N0z*dz;
   double delz = Nz*dz;

   printf("Domain (Nr, Nz): (%d, %d)\n r0: %f\n delr: %f\n z0: %f\n  delz: %f\n", Nr,Nz, r0,delr,z0,delz);
//------------------------------------------------------------------------------

  srand(theDomain->rank);
  rand();
  
  double r, z, phi;
  struct tracerList *theList = theDomain->theTracers;
  struct tracer *tr = theList->head;
  while( tr != NULL){
    r = getRandIn(r0, delr);
    z = getRandIn(z0, delz);
    phi = getRandIn(0, phi_max);
    //printf("Got this far\n");
    tr->R = r; tr->Z = z; tr->Phi = phi;
    tr->Type = 0; tr->rmFlag = 0;  
/*  if( tr->next==NULL )
       printf("Tr->Next does NOT exist\n");
    else
       printf("Tr->Next DOES exist\n");
*/
    tr = tr->next;
  }
  printTracerCoords( theDomain );
}

/*
void distributeTracers( struct domain *theDomain){
   //Distributes tracers to their appropriate processes
  //(rather tells processors to only track its tracers)

   int Num_R = theDomain->theParList.Num_R;
   int Num_Z = theDomain->theParList.Num_Z;
   int * dim_rank = theDomain->dim_rank;
   int * dim_size = theDomain->dim_size;
   double phi_max = theDomain->phi_max;

   int N0r = getN0( dim_rank[0]   , dim_size[0] , Num_R );
   int N1r = getN0( dim_rank[0]+1 , dim_size[0] , Num_R );

   int N0z = getN0( dim_rank[1]   , dim_size[1] , Num_Z );
   int N1z = getN0( dim_rank[1]+1 , dim_size[1] , Num_Z );

   double dr = 1./(double)Num_R;
   double rm = (double)N0r * dr;
   double rp = (double)N1r * dr;

   double dz = 1./(double)Num_Z;
   double zm = (double)N0z * dz;
   double zp = (double)N1z * dz;

   double xp[3] = {rp, phi_max, zp};
   double xm[3] = {rm, 0.0, zm};

   //Finish giving each process it's list of tracers
   
}
*/


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
   tr->Type  = type;
   if( c != NULL ){
     test_cell_vel( tr, c );
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
