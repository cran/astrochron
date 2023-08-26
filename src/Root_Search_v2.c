/* This code is a component of astrochron: An R Package for Astrochronology
 Copyright (C) 2018 Huaran Liu

 Contact Huaran Liu (huaranliu@gmail.com) for information on
 updates. 

 Legal Notice: Root_Search is a function for numerically searching for the roots for nonlinear equation from the boundary condition of Schinck and Guinasso 1975 model. Bisection method is  used here and the root is considered to be determined when the binary search interval is less than eps(defined in the program). 

 This program is distributed in the hope that it will be useful, 
 but WITHOUT ANY WARRANTY; without even the implied warranty of 
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE */


#include <stdint.h>
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <stdio.h>
#include <math.h>

SEXP  Root_Search(SEXP Gval)
{
 /* Input parameter*/
 double *G = REAL(Gval);
 
 /* Other parameters*/
 int i, N = 300;                     // Number of alphas 
 double Il, Iu, Im;                  // lower, upper, middle point
 /* double fl, fu, fm, rt =  0.0001; modified for CRAN compliance, SRM Aug. 24, 2023 */
 double fl, fm, rt =  0.0001;
 double eps = 10E-13;

 /* Output alphas */
 SEXP alpha = PROTECT(allocVector( REALSXP, N));
 double *a = REAL(alpha);

  for( i = 0; i< N; i++){
       Il = i*PI+eps;
       Iu = rt + PI;
       Im = (Iu+Il)/2;
       do{

          fl = Il/(*G)*cos(Il)/sin(Il)-Il*Il+1.0/4.0/(*G)/(*G);
 /*       fu = Iu/(*G)*cos(Iu)/sin(Iu)-Iu*Iu+1.0/4.0/(*G)/(*G); removed for CRAN compliance, SRM Aug. 24, 2023 */
          fm = Im/(*G)*cos(Im)/sin(Im)-Im*Im+1.0/4.0/(*G)/(*G);
          if( fl*fm > 0.0){
              Il = Im;}
          else{Iu = Im;}
          Im = (Iu+Il)/2;
          rt = Im;
       }while(fabs(Iu-Il)>= eps);
       a[i] = rt;
 }
      UNPROTECT(1);
      return alpha;
}
