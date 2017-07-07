
#include "../paul.h"

static double gamma_law = 0.0;

void setICparams( struct domain * theDomain ){
   gamma_law = theDomain->theParList.Adiabatic_Index;
}

void initial( double * prim , double * x ){
   
   double a = 0.2;
   double x0 = -0.5;
   
   double r   = x[0];
   double phi = x[1];
   double z   = x[2];
   double xx = r*cos(phi);

   double sig = 0.4;
   double frac = ( xx-x0 )/sig;
   if( -sig<(xx-x0) && (xx-x0)<sig ){ frac = 1.0; }

   prim[RHO] = 1.0 + a*( 1-pow(frac,4.0) );;
   prim[PPP] = pow(prim[RHO],gamma_law);
   prim[URR] = 0.0;
   prim[UPP] = 0.0;
   prim[UZZ] = 0.0;
   if( NUM_N>0 ) prim[NUM_C] = 0.0;
}
