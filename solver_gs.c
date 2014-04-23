///////////////////////////////////////////////////////////////////////////////
///
/// \file   solver_gs.c
///
/// \brief  Gauss Seidel equation sovler.
///
/// \author Mingang Jin, Qingyan Chen
///         Purdue University
///         Jin56@purdue.edu, YanChen@purdue.edu
///         Wangda Zuo
///         University of Miami
///         W.Zuo@miami.edu
///
/// \date   04/02/2013
///
/// This file provides functions for solving the discretize equation by
/// iterations. In FFD, the diffusion equation and pressure equaion are solved
/// through iterative approaches.
///
///////////////////////////////////////////////////////////////////////////////
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "data_structure.h"
#include "solver_gs.h"
#include "boundary.h"
#include "utility.h"

FILE *file1;

///////////////////////////////////////////////////////////////////////////////
///\brief Gauss-Seidel solver for pressure
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param Type Type of variable
///\param x Pointer to variable
///
///\return No return needed
///////////////////////////////////////////////////////////////////////////////
void GS_P(PARA_DATA *para, REAL **var, int Type, REAL *x) {
  REAL *as = var[AS], *aw = var[AW], *ae = var[AE], *an = var[AN];
  REAL *ap = var[AP], *af = var[AF], *ab = var[AB], *b = var[B];
  int imax = para->geom->imax, jmax= para->geom->jmax;
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);  
  int i, j, k, it=0;
  float tmp1,tmp2,residual=1;
  REAL *flagp = var[FLAGP],*flagu = var[FLAGU],*flagv = var[FLAGV];
  REAL *flagw = var[FLAGW];

  while(residual>0.001 && it <5) {
    tmp1=0;
    tmp2=0.0000000001;
    it +=1;
    residual=0;
    
    for(i=1; i<=imax; i++)
      for(j=1; j<=jmax; j++)
        for(k=1; k<=kmax; k++) {
          
          if (flagp[FIX(i,j,k)]==1) continue;

          x[FIX(i,j,k)] = ( ae[FIX(i,j,k)]*x[FIX(i+1,j,k)] 
                          + aw[FIX(i,j,k)]*x[FIX(i-1,j,k)]
                          + an[FIX(i,j,k)]*x[FIX(i,j+1,k)]
                          + as[FIX(i,j,k)]*x[FIX(i,j-1,k)]
                          + af[FIX(i,j,k)]*x[FIX(i,j,k+1)]
                          + ab[FIX(i,j,k)]*x[FIX(i,j,k-1)]
                          + b[FIX(i,j,k)] ) / ap[FIX(i,j,k)];
        }

    for(i=imax; i>=1; i--)
      for(j=1; j<=jmax; j++)
        for(k=1; k<=kmax; k++) {
          
          if (flagp[FIX(i,j,k)]==1) continue;

          x[FIX(i,j,k)] = (  ae[FIX(i,j,k)]*x[FIX(i+1,j,k)] 
                          + aw[FIX(i,j,k)]*x[FIX(i-1,j,k)]
                          + an[FIX(i,j,k)]*x[FIX(i,j+1,k)]
                          + as[FIX(i,j,k)]*x[FIX(i,j-1,k)]
                          + af[FIX(i,j,k)]*x[FIX(i,j,k+1)]
                          + ab[FIX(i,j,k)]*x[FIX(i,j,k-1)]
                          + b[FIX(i,j,k)] ) / ap[FIX(i,j,k)];
        }


    for(j=1; j<=jmax; j++)
      for(i=1; i<=imax; i++)
        for(k=1; k<=kmax; k++) {
          
          if (flagp[FIX(i,j,k)]==1) continue;

          x[FIX(i,j,k)] = (  ae[FIX(i,j,k)]*x[FIX(i+1,j,k)] 
                          + aw[FIX(i,j,k)]*x[FIX(i-1,j,k)]
                          + an[FIX(i,j,k)]*x[FIX(i,j+1,k)]
                          + as[FIX(i,j,k)]*x[FIX(i,j-1,k)]
                          + af[FIX(i,j,k)]*x[FIX(i,j,k+1)]
                          + ab[FIX(i,j,k)]*x[FIX(i,j,k-1)]
                          + b[FIX(i,j,k)] ) / ap[FIX(i,j,k)];
        }

     for(j=jmax; j>=1; j--)
       for(i=1; i<=imax; i++)
         for(k=1; k<=kmax; k++) {
           if (flagp[FIX(i,j,k)]==1) continue;
           
           x[FIX(i,j,k)] = (  ae[FIX(i,j,k)]*x[FIX(i+1,j,k)] 
                          + aw[FIX(i,j,k)]*x[FIX(i-1,j,k)]
                          + an[FIX(i,j,k)]*x[FIX(i,j+1,k)]
                          + as[FIX(i,j,k)]*x[FIX(i,j-1,k)]
                          + af[FIX(i,j,k)]*x[FIX(i,j,k+1)]
                          + ab[FIX(i,j,k)]*x[FIX(i,j,k-1)]
                          + b[FIX(i,j,k)] ) / ap[FIX(i,j,k)];
         }
 

    FOR_EACH_CELL
      if (flagp[FIX(i,j,k)]==1) continue;
      tmp1 += (float) fabs(ap[FIX(i,j,k)]*x[FIX(i,j,k)] 
          - ae[FIX(i,j,k)]*x[FIX(i+1,j,k)] - aw[FIX(i,j,k)]*x[FIX(i-1,j,k)]
          - an[FIX(i,j,k)]*x[FIX(i,j+1,k)] - as[FIX(i,j,k)]*x[FIX(i,j-1,k)]
          - af[FIX(i,j,k)]*x[FIX(i,j,k+1)] - ab[FIX(i,j,k)]*x[FIX(i,j,k-1)]
          - b[FIX(i,j,k)]);
          
      tmp2 += (float) fabs(ap[FIX(i,j,k)]*x[FIX(i,j,k)]);
    END_FOR
  
    residual  = (tmp1 /tmp2);
   
    // para->prob->resp=residual;
 
    // para->prob->iterp +=1; 
  }

} // End of GS_P( )




///////////////////////////////////////////////////////////////////////////////
///\brief Gauss-Seidel solver
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param flag Pointer to the cell property flag
///\param x Pointer to variable
///
///\return No return needed
///////////////////////////////////////////////////////////////////////////////
void Gauss_Seidel(PARA_DATA *para, REAL **var, REAL *flag, REAL *x) {
  REAL *as = var[AS], *aw = var[AW], *ae = var[AE], *an = var[AN];
  REAL *ap = var[AP], *af = var[AF], *ab = var[AB], *b = var[B];
  int imax = para->geom->imax, jmax= para->geom->jmax;
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);  
  int i, j, k, it=0;
  float tmp1,tmp2,residual=1;

  while(residual>0.001 && it<5) {
    tmp1=0;
    tmp2=0.0000000001;
    it +=1;
    residual=0;

    for(i=1; i<=imax; i++)
      for(j=1; j<=jmax; j++)
        for(k=1; k<=kmax; k++) {
          if (flag[FIX(i,j,k)]>=0) continue;
          
          x[FIX(i,j,k)] = (  ae[FIX(i,j,k)]*x[FIX(i+1,j,k)] 
                            + aw[FIX(i,j,k)]*x[FIX(i-1,j,k)]
                            + an[FIX(i,j,k)]*x[FIX(i,j+1,k)]
                            + as[FIX(i,j,k)]*x[FIX(i,j-1,k)]
                            + af[FIX(i,j,k)]*x[FIX(i,j,k+1)]
                            + ab[FIX(i,j,k)]*x[FIX(i,j,k-1)]
                            + b[FIX(i,j,k)] ) / ap[FIX(i,j,k)];
        }

    
    for(i=imax; i>=1; i--)
      for(j=1; j<=jmax; j++)
        for(k=1; k<=kmax; k++) {
          if (flag[FIX(i,j,k)]>=0) continue;
          
          x[FIX(i,j,k)] = (  ae[FIX(i,j,k)]*x[FIX(i+1,j,k)] 
                            + aw[FIX(i,j,k)]*x[FIX(i-1,j,k)]
                            + an[FIX(i,j,k)]*x[FIX(i,j+1,k)]
                            + as[FIX(i,j,k)]*x[FIX(i,j-1,k)]
                            + af[FIX(i,j,k)]*x[FIX(i,j,k+1)]
                            + ab[FIX(i,j,k)]*x[FIX(i,j,k-1)]
                            + b[FIX(i,j,k)] ) / ap[FIX(i,j,k)];
        }



    for(j=1; j<=jmax; j++)
      for(i=1; i<=imax; i++)
        for(k=1; k<=kmax; k++) {
          if (flag[FIX(i,j,k)]>=0) continue;
          
          x[FIX(i,j,k)] = (  ae[FIX(i,j,k)]*x[FIX(i+1,j,k)] 
                            + aw[FIX(i,j,k)]*x[FIX(i-1,j,k)]
                            + an[FIX(i,j,k)]*x[FIX(i,j+1,k)]
                            + as[FIX(i,j,k)]*x[FIX(i,j-1,k)]
                            + af[FIX(i,j,k)]*x[FIX(i,j,k+1)]
                            + ab[FIX(i,j,k)]*x[FIX(i,j,k-1)]
                            + b[FIX(i,j,k)] ) / ap[FIX(i,j,k)];
        }

    for(j=jmax; j>=1; j--)
      for(i=1; i<=imax; i++)
        for(k=1; k<=kmax; k++) {
          
          if (flag[FIX(i,j,k)]>=0) continue;
          
          x[FIX(i,j,k)] = (  ae[FIX(i,j,k)]*x[FIX(i+1,j,k)] 
                            + aw[FIX(i,j,k)]*x[FIX(i-1,j,k)]
                            + an[FIX(i,j,k)]*x[FIX(i,j+1,k)]
                            + as[FIX(i,j,k)]*x[FIX(i,j-1,k)]
                            + af[FIX(i,j,k)]*x[FIX(i,j,k+1)]
                            + ab[FIX(i,j,k)]*x[FIX(i,j,k-1)]
                            + b[FIX(i,j,k)] ) / ap[FIX(i,j,k)];
        }

     
    FOR_EACH_CELL
      if (flag[FIX(i,j,k)]>=0) continue;
      tmp1 += (float) fabs(ap[FIX(i,j,k)]*x[FIX(i,j,k)] 
          - ae[FIX(i,j,k)]*x[FIX(i+1,j,k)] - aw[FIX(i,j,k)]*x[FIX(i-1,j,k)]
          - an[FIX(i,j,k)]*x[FIX(i,j+1,k)] - as[FIX(i,j,k)]*x[FIX(i,j-1,k)]
          - af[FIX(i,j,k)]*x[FIX(i,j,k+1)] - ab[FIX(i,j,k)]*x[FIX(i,j,k-1)]
          - b[FIX(i,j,k)]);
      tmp2 += (float) fabs(ap[FIX(i,j,k)]*x[FIX(i,j,k)]);

    END_FOR
    
    residual  = (tmp1 /tmp2);
  }
   

}// End of Gauss-Seidel( )



