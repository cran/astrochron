/* This code is a component of astrochron: An R Package for Astrochronology
 Copyright (C) 2018 Huaran Liu

 Contact Huaran Liu (huaranliu@gmail.com) for information on
 updates. 

 Legal Notice: Impulse_Response is a function for calculating the impulse response using 1-D 
 advection diffusion model with Neymann boundary condition. For equations and derivations
 see Schink and Guinasso 1975 

 This program is distributed in the hope that it will be useful, 
 but WITHOUT ANY WARRANTY; without even the implied warranty of 
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE */




#include <stdint.h>
#include <math.h>
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <Rdefines.h>

SEXP Impulse_Response(SEXP delta_t, SEXP NT, SEXP Gval, SEXP alpha, SEXP verbose){

        /* transform input paramters and array from SEXP to C data type */
        double *dt = REAL(delta_t);
        int  *nt = INTEGER(NT);
        double *G = REAL(Gval);
	PROTECT(alpha = AS_NUMERIC(alpha));
        double *a; a  = NUMERIC_POINTER(alpha);
        int *verb = INTEGER(verbose);
        
        /* parameters needed for calculation */
        int  N = length(alpha);                           // Number of alphas
        int nx = (int)1.0/(*dt);                          // length of mix layer depth
        int depth_nt = 1000;                              // length of total record/output
        double xx[nx], tt[depth_nt];                      
	int i, j, k;
	double fsum, fsin, fcos, fexp,beta[N];
     
        /* output vector fc */
        SEXP fc = PROTECT(allocVector(REALSXP,depth_nt));
        double *f = REAL(fc);

        /* Initialize of vector fc to zeros */
        for(i = 0; i < depth_nt; i++){
            REAL(fc)[i] = 0.0;
        }
        
        /* print out experiment parameters */
/*      if(*verb){
        printf("\n verbose = %d", *verb);
        printf("\n -------   Impulse Response Function Calculation   ------- \n");
        printf("\n Unit in depth xx: 1 (Time it takes to travel through mix layer)");
        printf("\n Unit in time tt: 1 (Mix layer depth)");
        printf("\n dx = dt = %3.2f", *dt);
        printf("\n Time steps after deposition : nt = %d", *nt);
        printf("\n Length of total record / output series: depth_nt = %d", depth_nt);
        printf("\n Number of alpha calculated from boundary conditions(2N) : N = %d \n", N);
        } */

        /* Start Time tt Loop ... */
        for( j = 0; j < *nt; j++){
            tt[j] = ( (double)j+1.0)*(*dt);

             /* Start Depth xx Loop ... */ 

               /* Loop through outside mix layer (Replacement) */
               for( i = depth_nt-1; i >= nx; i--){ 
                f[i] = f[i-1];}

               /* Loop through mix layer */
               for( i = 0; i < nx; i++){
                   xx[i] = ((double)i+1.0)*(*dt);

                   /* Loop through all alphas */
                   fsum = 0.0;
	           for (k = 0; k < N; k++){ 

                      beta[k] = 0.5+0.5/(*G)/pow(a[k],2)+0.125/pow(*G,2)/pow(a[k],2);
                      fsin = 0.5/(*G)/a[k]*sin( a[k]*xx[i]);
                      fcos = cos( a[k]*xx[i]);
                      fexp = exp(-*G*pow(a[k],2)*tt[j]);
                      fsum = fsum + 1./beta[k]*(fcos + fsin )*fexp;

                      fsin = -0.5/(*G)/a[k]*sin( -a[k]*xx[i]);
                      fcos = cos(-a[k]*xx[i]);
                      fexp = exp(-*G*pow(-a[k],2)*tt[j]);
                      fsum = fsum + 1./beta[k]*(fcos + fsin )*fexp;
                   }
                   REAL(fc)[i] = exp(xx[i]/2.0/(*G)-tt[j]/4.0/(*G))*fsum;
	       }
	}

        UNPROTECT(2);        
        return(fc);

}
