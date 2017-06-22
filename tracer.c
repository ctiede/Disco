#include "paul.h"

void setTracerParams( struct domain * theDomain){

   theDomain->Ntr = 500;
}

void initializeTracers( struct domain *theDomain ){

   struct param_list theParamList = theDomain->theParList;
   double rmin = theParamList.rmin;
   double rmax = theParamList.rmax;
   double dr   = rmax - rmin;
   double zmin = theParamList.zmin;
   double zmax = theParamList.zmax;
   double dz   = zmax - zmin;
   double phimax = theParamList.phimax;

   srand(0);
   int Ntr = theDomain->Ntr;
   int i;
   for( i=0; i<Ntr; ++i ){
	struct tracer *tr = theDomain->theTracers+i;
        double r = rmin + ((double)rand()/(double)RAND_MAX)*dr;
	double z = zmin + ((double)rand()/(double)RAND_MAX)*dz;
	double phi = ((double)rand()/(double)RAND_MAX)*phimax;
	tr->R = r; tr->Z = z; tr->Phi = phi;
   }

}


int check_phi(double phi, double phip, double phim, double phi_max){

   if( phip > phim ){
	if( phim < phi && phi < phip){ return 1; }
   }
   else{
	if( phim < phi && phi < phi_max){
		if( 0.0 < phi && phi < phip){ return 1; }
	}
   }
   return 0;

}

int check_in_cell(struct tracer *Tr, double *xp, double *xm, double phi_max){

   double r   = Tr->R;
   double phi = Tr->Phi;
   double z   = Tr->Z;

   if( xm[0] < r && r < xp[0] ){
	if( xm[2] < z && z < xp[2]){
	     if( check_phi(phi, xp[1], xm[1], phi_max) ){
		return 1;
	     }
	}
   }
   return 0;

}

void get_local_vel(struct tracer *aTracer, struct cell *c){

   double vr, om, vz; 

   if( c != NULL ){
      vr = c->prim[URR];
      om = c->prim[UPP];
      vz = c->prim[UZZ];
   } else{
      vr = 0;
      om = 0;
      vz = 0;
   }
   aTracer->Vr    = vr;
   aTracer->Omega = om;
   aTracer->Vz    = vz;

}

struct cell * get_tracer_cell(struct domain *theDomain, struct tracer *aTracer){

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
		if( check_in_cell(aTracer, xp, xm, phi_max) ){
			return c;
		}
	   }
	}
   }
   return NULL;

}

void clean_pi(struct domain * );

void moveTracers(struct domain *theDomain, struct tracer *aTracer, double dt){

   double rmin = theDomain->theParList.rmin;
   double rmax = theDomain->theParList.rmax;
   double zmin = theDomain->theParList.zmin;
   double zmax = theDomain->theParList.zmax;

   double r = aTracer->R;
   double phi = aTracer->Phi;
   double z = aTracer->Z;

   r += aTracer->Vr*dt;
   if( r > rmax ) r = rmax;
   if( r < rmin ) r = rmin;

   z += aTracer->Vz*dt;
   if( z > zmax ) z = zmax;
   if( z < zmin ) z = zmin; 

   aTracer->R = r;
   aTracer->Z = z;

   phi += aTracer->Omega*dt;
   aTracer->Phi = phi; 
   clean_pi(theDomain);

}

void updateTracers(struct domain *theDomain, double dt){
 
   int Ntr = theDomain->Ntr;
   int tr;
   for( tr=0 ; tr<Ntr ; ++tr){
	struct cell *c = get_tracer_cell( theDomain, theDomain->theTracers + tr);
	get_local_vel( theDomain->theTracers + tr, c);
	moveTracers(theDomain, theDomain->theTracers + tr, dt);
   }
   
}
