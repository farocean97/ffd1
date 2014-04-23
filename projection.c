///////////////////////////////////////////////////////////////////////////////
///
/// \file   projection.c
///
/// \brief  Subroutines for performing pressure projection used in FFD
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
/// This file provides functions for performing pressure projection in FFD, which
/// is realized in two steps,including solving the pressure equation and performing
/// the velocity correction.
///
///////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include "math.h"
#include "data_structure.h"
#include "boundary.h"
#include "projection.h"
#include "solver_gs.h"
#include "utility.h"

FILE *file1;

///////////////////////////////////////////////////////////////////////////////
/// \brief Project the velocity
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param BINDEX Pointer to boundary index
///
///\return No return needed
///////////////////////////////////////////////////////////////////////////////
void projection(PARA_DATA *para, REAL **var, int **BINDEX) {
  int i, j, k;
  int imax = para->geom->imax, jmax = para->geom->jmax;
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  REAL dt= para->mytime->dt;
  REAL *x = var[X], *y = var[Y], *z = var[Z];
  REAL *gx = var[GX], *gy = var[GY], *gz = var[GZ]; 
  REAL *u = var[VX], *v = var[VY], *w = var[VZ];  
  REAL *p = var[IP], *b = var[B], *ap = var[AP], *ab = var[AB], *af = var[AF];
  REAL *ae = var[AE], *aw =var[AW], *an = var[AN], *as = var[AS];
  REAL *pp = var[PP];
  REAL dxe,dxw, dyn,dys,dzf,dzb,Dx,Dy,Dz;
  REAL zv=para->geom->zv;
  REAL residual = 1.0;  
  REAL *flagp = var[FLAGP],*flagu = var[FLAGU],*flagv = var[FLAGV],*flagw = var[FLAGW];


 // Specify the equation coefficients
  FOR_EACH_CELL

    dxe = x[FIX(i+1,j,k)]-x[FIX(i,j,k)];
    dxw = x[FIX(i,j,k)]-x[FIX(i-1,j,k)];
    dyn = y[FIX(i,j+1,k)]-y[FIX(i,j,k)];
    dys = y[FIX(i,j ,k)]-y[FIX(i,j-1,k)];
    dzf = z[FIX(i,j ,k+1)]-z[FIX(i,j,k)];
    dzb = z[FIX(i,j ,k)]-z[FIX(i,j,k-1)];
    Dx = gx[FIX(i,j,k)]-gx[FIX(i-1,j,k)];
    Dy = gy[FIX(i,j ,k)]-gy[FIX(i,j-1,k)];
    Dz = gz[FIX(i,j ,k)]-gz[FIX(i,j,k-1)];
 
    ae[FIX(i,j,k)] = Dy*Dz/dxe;
    aw[FIX(i,j,k)] = Dy*Dz/dxw;      
    an[FIX(i,j,k)] = Dx*Dz/dyn;
    as[FIX(i,j,k)] = Dx*Dz/dys;
    af[FIX(i,j,k)] = Dx*Dy/dzf;
    ab[FIX(i,j,k)] = Dx*Dy/dzb;

    b[FIX(i,j,k)] =  Dy*Dz/dt*u[FIX(i-1,j,k)]-Dy*Dz/dt*u[FIX(i,j,k)] 
                    + Dx*Dz/dt*v[FIX(i,j-1,k)]-Dx*Dz/dt*v[FIX(i,j,k)] 
                    + Dx*Dy/dt*w[FIX(i,j,k-1)]-Dx*Dy/dt*w[FIX(i,j,k)];

  END_FOR

  set_bnd_pressure(para, var, p,BINDEX); 

  FOR_EACH_CELL    
    ap[FIX(i,j,k)] = ae[FIX(i,j,k)] + aw[FIX(i,j,k)] + as[FIX(i,j,k)]
                   + an[FIX(i,j,k)] + af[FIX(i,j,k)] + ab[FIX(i,j,k)];
  END_FOR

  GS_P(para, var, IP, p);  

  vel_correction(para,var, BINDEX);

} // End of projecttion( )


///////////////////////////////////////////////////////////////////////////////
///\brief Correct the velocity based on pressure solution
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param BINDEX Pointer to boundary index
///
///\return No return needed
///////////////////////////////////////////////////////////////////////////////
void vel_correction(PARA_DATA *para, REAL **var,int **BINDEX) {
  int i, j, k;
  int imax = para->geom->imax, jmax = para->geom->jmax;
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  REAL dt= para->mytime->dt;
  REAL *x = var[X], *y = var[Y], *z = var[Z];
  REAL *gx = var[GX], *gy = var[GY], *gz = var[GZ];
  REAL *u = var[VX], *v = var[VY], *w = var[VZ]; 
  REAL *p = var[IP];
  REAL *flagp = var[FLAGP],*flagu = var[FLAGU];
  REAL *flagv = var[FLAGV],*flagw = var[FLAGW];
  REAL zv=para->geom->zv;

  set_bnd_pressure(para, var, p,BINDEX); 

  FOR_ALL_CELL
   
    if(flagu[FIX(i,j,k)]>=0) u[FIX(i,j,k)]=u[FIX(i,j,k)];
    else u[FIX(i,j,k)] = u[FIX(i,j,k)]- dt
                         *(p[FIX(i+1,j,k)]/(x[FIX(i+1,j,k)]-x[FIX(i,j,k)])
                          -p[FIX(i,j,k)]/(x[FIX(i+1,j,k)]-x[FIX(i,j,k)]));

    if(flagv[FIX(i,j,k)]>=0 ) v[FIX(i,j,k)]=v[FIX(i,j,k)];
    else {
          v[FIX(i,j,k)] = v[FIX(i,j,k)]- dt
                         *(p[FIX(i,j+1,k)]/ (y[FIX(i,j+1,k)]-y[FIX(i,j,k)])
                          -p[FIX(i,j,k)]/ (y[FIX(i,j+1,k)]-y[FIX(i,j,k)]));
    }

    if (flagw[FIX(i,j,k)]>=0) w[FIX(i,j,k)]=w[FIX(i,j,k)];
    else w[FIX(i,j,k)] = w[FIX(i,j,k)]- dt
                         *(p[FIX(i,j,k+1)] / (z[FIX(i,j,k+1)]-z[FIX(i,j,k)])
                          -p[FIX(i,j,k)] / (z[FIX(i,j,k+1)]-z[FIX(i,j,k)]));

  END_FOR

}
