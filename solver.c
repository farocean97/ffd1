///////////////////////////////////////////////////////////////////////////////
///
/// \file   solver.c
///
/// \brief  The solver of FFD
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
/// This file provides functions for performing the computation of solving the
/// governing equation.FFD solves the set of equations sequentially.First the
/// \c vel_step() is applied to solve the momentum equation, and then energy 
/// equation \c temp_step(), and species transport \c den_step().
///
///////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include "data_structure.h"
#include "solver.h"
#include "write_data.h"
#include "diffusion.h"
#include "projection.h"
#include "advection.h"
#include "timing.h"
#include "solver_gs.h"
#include "boundary.h"
#include "utility.h"
#include "chen_zero_equ_model.h"


///////////////////////////////////////////////////////////////////////////////
///\brief Entrance of the FFD solver
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param BINDEX Pointer to boundary index
///
///\return No return needed
///////////////////////////////////////////////////////////////////////////////
void FFD_solver(PARA_DATA *para, REAL **var,int **BINDEX) {
  int imax = para->geom->imax, jmax = para->geom->jmax;
  int kmax = para->geom->kmax;
  int size = (imax+2) * (jmax+2) * (kmax+2);
  int t_step = 0, t_output = para->mytime->t_output;
  REAL t_steady = para->mytime->t_steady;
  REAL dt = para->mytime->dt;
  REAL *u = var[VX], *v = var[VY], *w = var[VZ];
  REAL *den = var[DEN], *temp = var[TEMP];
  int tsize=3;
  int  IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  REAL *x = var[X], *y = var[Y], *z = var[Z];
  REAL *gx = var[GX], *gy = var[GY], *gz = var[GZ];
  REAL *vsx=var[VSX];
  //int cal_mean = para->outp->cal_mean;

  lengthscale(para,var); //zero equation length scale
  
  while( para->mytime->t_step < t_output) {
    vis_turb(para, var);      
    vel_step(para, var, BINDEX);
    temp_step(para, var,BINDEX);
    timing(para);

  }
  // write out the result
  write_unsteady(para, var, "unsteady");
  write_Tecplot(para, var, "result",BINDEX);
  para->prob->output = 1;

} // End of FFD_solver( ) 


///////////////////////////////////////////////////////////////////////////////
///\brief step for solving temperature
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param BINDEX Pointer to boundary index
///
///\return No return needed
///////////////////////////////////////////////////////////////////////////////
void temp_step(PARA_DATA *para, REAL **var,int **BINDEX){

  REAL *T = var[TEMP], *T0 = var[TMP1],*T1 = var[TMP2];            
  REAL *flagp = var[FLAGP];
  int j=para->geom->jplume+1;
  int k=para->geom->kplume;
  int imax=para->geom->imax, jmax=para->geom->jmax;
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);

  set_bnd_temp(para, var, TEMP, T,BINDEX); 

  advection(para, var, TEMP,flagp, T0, T, BINDEX);

  if(para->prob->plume_mod==1) plume_thermal(para,var,BINDEX);

  psi_conservation(para, var,T0,T,BINDEX);

  diffusion(para, var, TEMP, T, T0, BINDEX);

} // End of temp_step( )


///////////////////////////////////////////////////////////////////////////////
///\brief step for solving specices transport
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param BINDEX Pointer to boundary index
///
///\return No return needed
///////////////////////////////////////////////////////////////////////////////
void den_step(PARA_DATA *para, REAL **var,int **BINDEX){
  REAL *den = var[DEN], *den0 = var[TMP1];
  REAL *flagp = var[FLAGP];
      
  set_bnd_density(para, var,DEN,den,BINDEX);
  advection(para, var, DEN,flagp, den0, den,BINDEX); 
  psi_conservation(para, var,den0,den,BINDEX);
  diffusion(para, var, DEN, den, den0,BINDEX); 


} // End of den_step( )



///////////////////////////////////////////////////////////////////////////////
///\brief Step for solving the velocity
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param BINDEX Pointer to boundary index
///
///\return No return needed
///////////////////////////////////////////////////////////////////////////////
void vel_step(PARA_DATA *para, REAL **var,int **BINDEX) {
  REAL *u  = var[VX],  *v  = var[VY],    *w  = var[VZ];
  REAL *u0 = var[TMP1], *v0 = var[TMP2], *w0 = var[TMP3];
  REAL *u1 = var[VXM], *v1 = var[VYM], *w1 = var[VZM];
  REAL *p=var[IP];
  REAL *flagu = var[FLAGU],*flagv = var[FLAGV],*flagw = var[FLAGW];
  REAL *b=var[B];
  int i=para->geom->iplume;
  int j=para->geom->jplume;
  int k=para->geom->kplume;
  int imax=para->geom->imax, jmax=para->geom->jmax;
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
        

  set_bnd(para, var, VX, u, BINDEX);
  set_bnd(para, var, VY, v, BINDEX);
  set_bnd(para, var, VZ, w, BINDEX);


  advection(para, var, VX, flagu, u0, u, BINDEX);
  advection(para, var, VY, flagv, v0, v, BINDEX);       
  advection(para, var, VZ, flagw, w0, w, BINDEX); 


  diffusion(para, var, VX, u1, u0, BINDEX);   
  diffusion(para, var, VY, v1, v0, BINDEX); 
  diffusion(para, var, VZ, w1, w0, BINDEX); 

  if(para->bc->NBOUT!=0) mass_conservation(para, var,BINDEX);
  swapuvw(para,var);

  projection(para, var,BINDEX);

  if(para->prob->plume_mod==1) plume(para,var,BINDEX);

} // End of vel_step( )



///////////////////////////////////////////////////////////////////////////////
///\brief Entrance of iterative equation solver
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param var_type The type of variable
///\param psi Pointer to the solving variable
///
///\return No return needed
///////////////////////////////////////////////////////////////////////////////
void equ_solver(PARA_DATA *para, REAL **var, int var_type, REAL *psi) {
  REAL *flagp = var[FLAGP],*flagu = var[FLAGU],*flagv = var[FLAGV];
  REAL *flagw = var[FLAGW];

  switch (var_type) {
    case VX:
      Gauss_Seidel(para, var, flagu, psi);
    break;
    case VY:
      Gauss_Seidel(para, var, flagv, psi);
    break;
    case VZ:
      Gauss_Seidel(para, var, flagw, psi);
    break;
    case TEMP:
    case IP:
    case DEN:
      GS_P(para, var, TEMP, psi);
    break;
  }

}// End of equ_solver()