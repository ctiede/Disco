#include "paul.h"

void setTracerParams( struct domain * theDomain){

   theDomain->Ntr = 20;
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
        tr->Type = 0;
   }

}

void clean_pi_tr( struct domain *theDomain ){
   int Ntr = theDomain->Ntr;
   double phi_max = theDomain->phi_max;
   int i;
   for( i=0; i<Ntr; ++i){
      struct tracer *tr = theDomain->theTracers + i;
      double phi = tr->Phi;
      while( phi > phi_max ){ phi -= phi_max; }
      tr->Phi = phi;
   }

}


int check_phi(double phi, double phip, double phim, double phi_max){
/*
   if( phip > phim ){
	if( phim < phi && phi < phip){ return 1; }
   }
   else{
	if( phim<phi && phi<phi_max ){ return 1; }
	else if( 0.0<phi && phi<phip ){ return 1; }
   }
*/

  if( phim < phi && phi < phip ){ return 1; }
  else{
   return 0;
  }
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

void test_cell_vel( struct tracer *tr, struct cell *c ){

  int check1=0, check2=0, check3=0;
  if( isnan(c->prim[URR]) ){
    tr->Vr 	  = 0;
    tr->Omega = 0;
    tr->Vz    = 0;
    tr->Type  = 3;
    check1 = 1;
   }/*
   if( c->prim[UPP] != c->prim[UPP] ){
	tr->Vr    = 0;
        tr->Omega = 0;
        tr->Vz    = 0;
        tr->Type  = 4;
	check2 = 1;
   }
   if( c->prim[UZZ] != c->prim[UZZ] ){
	tr->Vr    = 0;
        tr->Omega = 0;
        tr->Vz    = 0;
        tr->Type  = 5;
	check3 = 1;
   }

  if( check1==1 && check2==1 && check3==1){ tr->Type = 6; }
*/
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

//void clean_pi(struct domain * );

void moveTracers(struct domain *theDomain, struct tracer *aTracer, double dt){

   double rmin = theDomain->theParList.rmin;
   double rmax = theDomain->theParList.rmax;
   double zmin = theDomain->theParList.zmin;
   double zmax = theDomain->theParList.zmax;

   double r = aTracer->R;
   double phi = aTracer->Phi;
   double z = aTracer->Z;

  // double k1, k2;

   r += aTracer->Vr*dt;
   if( r > rmax ) r = 0.0/0.0;
   if( r < rmin ) r = 0.0/0.0;

   z += aTracer->Vz*dt;
   if( z > zmax ) z = 0.0/0.0;
   if( z < zmin ) z = 0.0/0.0;

   aTracer->R = r;
   aTracer->Z = z;

   phi += aTracer->Omega*dt;
   aTracer->Phi = phi;
   clean_pi_tr(theDomain);

}

/*
void tracer_RK_copy( struct tracer * tr ){
   tr->RK_r     = tr->R;
   tr->RK_phi   = tr->Phi;
   tr->RK_z	= tr-Z;
   tr->RK_vr	= tr->Vr;
   tr->RK_omega = tr->Omega;
   tr->RK_vz	= tr->Vz;
}

void tracer_RK_adjust( struct planet * pl , double RK ){
   tr->R     = (1.-RK)*tr->R     + RK*tr->RK_r;
   tr->Phi   = (1.-RK)*tr->Phi   + RK*tr->RK_phi;
   tr->Z     = (1.-RK)*tr->Z     + RK*tr->RK_z;
   tr->Vr    = (1.-RK)*tr->Vr    + RK*tr->RK_vr;
   tr->Omega = (1.-RK)*tr->Omega + RK*tr->RK_omega;
   tr->Vz    = (1.-RK)*tr->Vz    + RK*tr->RK_vz;
}
*/



void updateTracers(struct domain *theDomain, double dt){

   int Ntr = theDomain->Ntr;
   int tr;
   for( tr=0 ; tr<Ntr ; ++tr){
	struct cell *c = get_tracer_cell( theDomain, theDomain->theTracers + tr);
	get_local_vel( theDomain->theTracers + tr, c);
	//test_cell_vel( c,  );
        moveTracers(theDomain, theDomain->theTracers + tr, dt);
   }

}
