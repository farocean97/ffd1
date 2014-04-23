///////////////////////////////////////////////////////////////////////////////
///
/// \file   chen_zero_equ_model.c
///
/// \brief  Zero eqution model
///
/// \author Mingang Jin, Qingyan Chen
///         Purdue University
///         Jin56@purdue.edu, YanChen@purdue.edu
///         Wangda Zuo
///         University of Miami
///         W.Zuo@miami.edu
///
/// \date   04/02/2014
///
/// This file provides functions that used for calcuating the the turbulence
/// viscosity with zero equation model. The characteristic length in the model 
/// will be computed first with \c lengthscale(), and then the the turbulence
/// viscosity will be calculated.
///
///////////////////////////////////////////////////////////////////////////////
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>

#include "data_structure.h"

///////////////////////////////////////////////////////////////////////////////
/// \brief Computes lengthscale for zero equation model.
///
/// The distance from a fluid cell to the nearest boudnary cell will be
/// calculated.
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///
///\return void No return needed
///////////////////////////////////////////////////////////////////////////////
void lengthscale(PARA_DATA *para, REAL **var){
  int i,j,k;
  int p,q,r;
  REAL *x = var[X], *y = var[Y], *z = var[Z];
  REAL *gx = var[GX], *gy = var[GY], *gz = var[GZ];
  int imax = para->geom->imax, jmax = para->geom->jmax;
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  REAL *flagp=var[FLAGP],*flagu=var[FLAGU],*flagv=var[FLAGV],*flagw=var[FLAGW];
  REAL dx,dy,dz,dist=0;
  REAL  distmax=0;

  distmax=distmax> para->geom->Lx ? distmax : para->geom->Lx; 
  distmax=distmax> para->geom->Ly ? distmax : para->geom->Ly; 
  distmax=distmax> para->geom->Lz ? distmax : para->geom->Lz; 

  FOR_EACH_CELL
    if(flagp[FIX(i,j,k)]==1) continue;
    dist=distmax;
    
    //distant in x direction
    for(p=0;p<=imax+1;p++) {
      if(flagu[FIX(p,j,k)]>0) {
        dx=fabs(gx[FIX(p,j,k)]-x[FIX(i,j,k)]);
        dist = dist< dx ? dist : dx; 
      }
    }
     
    //distant in y direction   
    for(q=0;q<=jmax+1;q++)  {
      if(flagv[FIX(i,q,k)]>0) {
        dy=fabs(gy[FIX(i,q,k)]-y[FIX(i,j,k)]);
        dist = dist< dy ? dist : dy;
      }
    }
     
    //distant in y direction
    for(r=0;r<=kmax+1;r++) {
      if(flagw[FIX(i,j,r)]>0) {
        dz=fabs(gz[FIX(i,j,r)]-z[FIX(i,j,k)]);
        dist = dist< dz ? dist : dz;
      }
    }

    var[DIST][FIX(i,j,k)] = dist;

  END_FOR

} // End of lengthscale()


///////////////////////////////////////////////////////////////////////////////
/// \brief Formulation of zero equation model.
///
/// The model of zero equation model with coefficeint Chen_a=0.03847
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///
///\return REAL turbulenc viscosity
///////////////////////////////////////////////////////////////////////////////
REAL nu_t_chen_zero_equ(REAL u,REAL v,REAL w, REAL dist) {

  REAL chen_a=0.03847f;//0.03847f;
  REAL nu_t;
  nu_t = chen_a*dist*sqrt(u*u+v*v+w*w);
  return nu_t;
} // End of nu_t_chen_zero_eqn()


///////////////////////////////////////////////////////////////////////////////
/// \brief Computes tubulence viscosity for a fluid cell
///
/// The velocity magnitude at this cell will be calcuated before computing the 
/// turbulence visocity
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///
///\return void No return needed
///////////////////////////////////////////////////////////////////////////////
void vis_turb(PARA_DATA *para, REAL **var) {

  int i,j,k;
  REAL *u = var[VX], *v = var[VY], *w = var[VZ];
  REAL *vt=var[VT],*dist=var[DIST];
  int imax = para->geom->imax, jmax = para->geom->jmax;
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  REAL u0,v0,w0;
  
  FOR_EACH_CELL
    u0 = 0.5f*(u[FIX(i,j,k)]+u[FIX(i-1,j,k)]);
    v0 = 0.5f*(v[FIX(i,j,k)]+v[FIX(i,j-1,k)]);
    w0 = 0.5f*(w[FIX(i,j,k)]+w[FIX(i,j,k-1)]);
    vt[FIX(i,j,k)] = nu_t_chen_zero_equ(u0,v0,w0,dist[FIX(i,j,k)]);
  END_FOR

} // End of vis_turb()