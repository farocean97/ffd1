///////////////////////////////////////////////////////////////////////////////
///
/// \file   diffusion.c
///
/// \brief  Subroutines for solving the diffison equation
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
/// This file provides functions that used for solving the diffusion equation
/// in FFD. By calling the function \c diffusion(), the coefficients of the 
/// discretized equation can be specified through \c coef_diff(), and then the
/// equation solver will be applied to solve the equation.
///
///////////////////////////////////////////////////////////////////////////////
#include <stdio.h>
#include <stdlib.h>

#include "data_structure.h"
#include "diffusion.h"
#include "boundary.h"
#include "solver.h"
#include "utility.h"
#include "chen_zero_equ_model.h"

///////////////////////////////////////////////////////////////////////////////
///\brief Entrance of calculating diffusion equation
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param var_type Type of variable
///\param psi Pointer to the variable at current time step
///\param psi0 Pointer to the variable at previous time step
///\param BINDEX Pointer to boundary index
///
///\return void No return needed
///////////////////////////////////////////////////////////////////////////////
void diffusion(PARA_DATA *para, REAL **var, int var_type, 
               REAL *psi, REAL *psi0, int **BINDEX) {

  // Define the coefcients for euqations
  coef_diff(para, var, psi, psi0, var_type, BINDEX);

  // Solving the discretized equation
  equ_solver(para, var, var_type, psi);

  // Define the values of variables at B.C.
  set_bnd(para, var, var_type, psi, BINDEX);
          
} // End of diffusion( )


///////////////////////////////////////////////////////////////////////////////
///\brief Specifying the coefficients for equation
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param psi Pointer to the variable at current time step
///\param psi0 Pointer to the variable at previous time step
///\param var_type Type of variable
///\param BINDEX Pointer to boundary index
///
///\return void No return needed
///////////////////////////////////////////////////////////////////////////////
void coef_diff(PARA_DATA *para, REAL **var, REAL *psi, REAL *psi0, 
               int var_type, int **BINDEX) {
  int i, j, k;
  int imax = para->geom->imax, jmax = para->geom->jmax;
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  REAL *u=var[VX],*v=var[VY],*w=var[VZ];
  REAL *aw = var[AW], *ae = var[AE], *as = var[AS], *an = var[AN];
  REAL *af = var[AF], *ab = var[AB], *ap = var[AP], *ap0 = var[AP0], *b = var[B];
  REAL *x = var[X], *y = var[Y], *z = var[Z];
  REAL *gx = var[GX], *gy = var[GY], *gz = var[GZ];
  REAL *dist= var[DIST],*vt=var[VT];
  REAL *Temp = var[TEMP];
  REAL *d=var[DEN];
  REAL dxe, dxw, dyn, dys, dzf, dzb, Dx, Dy, Dz;
  REAL dt = para->mytime->dt, beta = para->prob->beta;
  REAL Temp_opt = para->prob->Temp_opt, gravx = para->prob->gravx;
  REAL gravy = para->prob->gravy, gravz = para->prob->gravz;
  REAL kapa;
  REAL *flagp = var[FLAGP],*flagu = var[FLAGU];
  REAL *flagv = var[FLAGV],*flagw = var[FLAGW];
  REAL Prt= para->prob->Prt;
  REAL dtemp=0.0f;
  REAL *s_plume=var[ER];



  switch (var_type)  {
    // coefficients for velocity u
    case VX: 
      //Specifying the value for viscosity.
      kapa = para->prob->nu; 

      FOR_U_CELL
        if(flagu[FIX(i,j,k)]>=0) continue;
        dxe = gx[FIX(i+1,j  ,k)] - gx[FIX(i  ,j,k)];
        dxw = gx[FIX(i  ,j  ,k)] - gx[FIX(i-1,j,k)];
        dyn =  y[FIX(i  ,j+1,k)] -  y[FIX(i  ,j,k)];
        dys=y[FIX(i,j,k)]-y[FIX(i,j-1,k)];
        dzf=z[FIX(i,j,k+1)]-z[FIX(i,j,k)];
        dzb=z[FIX(i,j,k)]-z[FIX(i,j,k-1)];
        Dx=x[FIX(i+1,j,k)]-x[FIX(i,j,k)];
        Dy=gy[FIX(i,j,k)]-gy[FIX(i,j-1,k)];
        Dz=gz[FIX(i,j,k)]-gz[FIX(i,j,k-1)];
        
        //zero equation model
        if(para->prob->tur_model == CHEN)  {
          aw[FIX(i,j,k)]= (kapa+vt[FIX(i,j,k)])*Dy*Dz/dxw;
          ae[FIX(i,j,k)]= (kapa+vt[FIX(i+1,j,k)])*Dy*Dz/dxe;
          an[FIX(i,j,k)]= (kapa+0.25f*(vt[FIX(i,j,k)]+vt[FIX(i+1,j,k)]
                          +vt[FIX(i,j+1,k)]+vt[FIX(i+1,j+1,k)]))*Dx*Dz/dyn;
          as[FIX(i,j,k)]= (kapa+0.25f*(vt[FIX(i,j,k)]+vt[FIX(i+1,j,k)]
                          +vt[FIX(i,j-1,k)]+vt[FIX(i+1,j-1,k)]))*Dx*Dz/dys;
          af[FIX(i,j,k)]= (kapa+0.25f*(vt[FIX(i,j,k)]+vt[FIX(i+1,j,k)]
                          +vt[FIX(i,j,k+1)]+vt[FIX(i+1,j,k+1)]))*Dx*Dy/dzf;
          ab[FIX(i,j,k)]= (kapa+0.25f*(vt[FIX(i,j,k)]+vt[FIX(i+1,j,k)]
                          +vt[FIX(i,j,k-1)]+vt[FIX(i+1,j,k-1)]))*Dx*Dy/dzb;
        }
        //laminar model
        else  { 
          aw[FIX(i,j,k)]= kapa*Dy*Dz/dxw;
          ae[FIX(i,j,k)]= kapa*Dy*Dz/dxe;
          an[FIX(i,j,k)]= kapa*Dx*Dz/dyn;
          as[FIX(i,j,k)]= kapa*Dx*Dz/dys;
          af[FIX(i,j,k)]= kapa*Dx*Dy/dzf;
          ab[FIX(i,j,k)]= kapa*Dx*Dy/dzb; 
        }
        
        ap0[FIX(i,j,k)]= Dx*Dy*Dz/dt;
        b[FIX(i,j,k)]= psi0[FIX(i,j,k)]*ap0[FIX(i,j,k)];
        //buoyancy forces && plume model correction
        if(para->prob->gravdir==GRAVX || para->prob->gravdir==GRAVXN ) {
          b[FIX(i,j,k)] -= beta*gravx *(0.5f*(Temp[FIX(i,j,k)]
                                       +Temp[FIX(i+1,j,k)])-Temp_opt)*Dx*Dy*Dz;
          if(para->prob->plume_mod==1) 
            b[FIX(i,j,k)] += s_plume[FIX(i,j,k)]/dt*Dx*Dy*Dz;//plume model  
        }

        ap[FIX(i,j,k)] = ap0[FIX(i,j,k)] + ae[FIX(i,j,k)] + aw[FIX(i,j,k)] 
                        + an[FIX(i,j,k)] + as[FIX(i,j,k)] + af[FIX(i,j,k)] 
                        + ab[FIX(i,j,k)];
      END_FOR

      set_bnd(para, var, var_type, psi, BINDEX);
      
      break;
    // coefficients for velocity v
    case VY:  
      kapa = para->prob->nu; 
      
    
      FOR_V_CELL
        if(flagv[FIX(i,j,k)]>=0) continue;
        dxe=x[FIX(i+1,j,k)]-x[FIX(i,j,k)];
        dxw=x[FIX(i,j,k)]-x[FIX(i-1,j,k)];
        dyn=gy[FIX(i,j+1,k)]-gy[FIX(i,j,k)];
        dys=gy[FIX(i,j,k)]-gy[FIX(i,j-1,k)];
        dzf=z[FIX(i,j,k+1)]-z[FIX(i,j,k)]; 
        dzb=z[FIX(i,j,k)]-z[FIX(i,j,k-1)];
        Dx=gx[FIX(i,j,k)]-gx[FIX(i-1,j,k)];
        Dy=y[FIX(i,j+1,k)]-y[FIX(i,j,k)];
        Dz=gz[FIX(i,j,k)]-gz[FIX(i,j,k-1)];

        if(para->prob->tur_model == CHEN)  {
          aw[FIX(i,j,k)]= (kapa+0.25f*(vt[FIX(i,j,k)]+vt[FIX(i,j+1,k)]
                          +vt[FIX(i-1,j,k)]+vt[FIX(i-1,j+1,k)]))*Dy*Dz/dxw;
          ae[FIX(i,j,k)]= (kapa+0.25f*(vt[FIX(i,j,k)]+vt[FIX(i,j+1,k)]
                          +vt[FIX(i+1,j,k)]+vt[FIX(i+1,j+1,k)]))*Dy*Dz/dxe;
          an[FIX(i,j,k)]= (kapa+vt[FIX(i,j+1,k)])*Dx*Dz/dyn;
          as[FIX(i,j,k)]= (kapa+vt[FIX(i,j,k)])*Dx*Dz/dys;
          af[FIX(i,j,k)]= (kapa+0.25f*(vt[FIX(i,j,k)]+vt[FIX(i,j+1,k)]
                          +vt[FIX(i,j,k+1)]+vt[FIX(i,j+1,k+1)]))*Dx*Dy/dzf;
          ab[FIX(i,j,k)]= (kapa+0.25f*(vt[FIX(i,j,k)]+vt[FIX(i,j+1,k)]
                          +vt[FIX(i,j,k-1)]+vt[FIX(i,j+1,k-1)]))*Dx*Dy/dzb;
        }
        else  {
          aw[FIX(i,j,k)]= kapa*Dy*Dz/dxw;
          ae[FIX(i,j,k)]= kapa*Dy*Dz/dxe;
          an[FIX(i,j,k)]= kapa*Dx*Dz/dyn;
          as[FIX(i,j,k)]= kapa*Dx*Dz/dys;
          af[FIX(i,j,k)]= kapa*Dx*Dy/dzf;
          ab[FIX(i,j,k)]= kapa*Dx*Dy/dzb; 
        }

        ap0[FIX(i,j,k)]= Dx*Dy*Dz/dt;
        b[FIX(i,j,k)]= psi0[FIX(i,j,k)]*ap0[FIX(i,j,k)];

        if(para->prob->gravdir==GRAVY || para->prob->gravdir==GRAVYN ) {

          b[FIX(i,j,k)] -= beta*gravy*(0.5f*(Temp[FIX(i,j,k)]
                                   +Temp[FIX(i,j+1,k)])-Temp_opt)*Dx*Dy*Dz;
          if(para->prob->plume_mod==1)
            b[FIX(i,j,k)] += s_plume[FIX(i,j,k)]/dt*Dx*Dy*Dz;//plume model 
        }
        ap[FIX(i,j,k)] = ap0[FIX(i,j,k)] + ae[FIX(i,j,k)] + aw[FIX(i,j,k)] 
                       + an[FIX(i,j,k)]  + as[FIX(i,j,k)] + af[FIX(i,j,k)] 
                       + ab[FIX(i,j,k)];
      END_FOR

      set_bnd(para, var, var_type, psi, BINDEX);

    break;
    // coefficients for velocity w
    case VZ: 
      kapa = para->prob->nu; 

      FOR_W_CELL
        if(flagw[FIX(i,j,k)]>=0) continue;
        dxe=x[FIX(i+1,j,k)]-x[FIX(i,j,k)];
        dxw=x[FIX(i,j,k)]-x[FIX(i-1,j,k)];
        dyn=y[FIX(i,j+1,k)]-y[FIX(i,j,k)];
        dys=y[FIX(i,j,k)]-y[FIX(i,j-1,k)];
        dzf=gz[FIX(i,j,k+1)]-gz[FIX(i,j,k)];
        dzb=gz[FIX(i,j,k)]-gz[FIX(i,j,k-1)];
        Dx=gx[FIX(i,j,k)]-gx[FIX(i-1,j,k)];
        Dy=gy[FIX(i,j,k)]-gy[FIX(i,j-1,k)];
        Dz=z[FIX(i,j,k+1)]-z[FIX(i,j,k)];

        if(para->prob->tur_model == CHEN) {
          aw[FIX(i,j,k)]= (kapa+0.25f*(vt[FIX(i,j,k)]+vt[FIX(i,j,k+1)]
                          +vt[FIX(i-1,j,k)]+vt[FIX(i-1,j,k+1)]))*Dy*Dz/dxw;
          ae[FIX(i,j,k)]= (kapa+0.25f*(vt[FIX(i,j,k)]+vt[FIX(i,j,k+1)]
                          +vt[FIX(i+1,j,k)]+vt[FIX(i+1,j,k+1)]))*Dy*Dz/dxe;
          an[FIX(i,j,k)]= (kapa+0.25f*(vt[FIX(i,j,k)]+vt[FIX(i,j,k+1)]
                          +vt[FIX(i,j+1,k)]+vt[FIX(i,j+1,k+1)]))*Dx*Dz/dyn;
          as[FIX(i,j,k)]= (kapa+0.25f*(vt[FIX(i,j,k)]+vt[FIX(i,j,k+1)]
                          +vt[FIX(i,j-1,k)]+vt[FIX(i,j-1,k+1)]))*Dx*Dz/dys;
          af[FIX(i,j,k)]= (kapa+vt[FIX(i,j,k+1)])*Dx*Dy/dzf;
          ab[FIX(i,j,k)]= (kapa+vt[FIX(i,j,k)])*Dx*Dy/dzb;
        }
        else {
          aw[FIX(i,j,k)]= kapa*Dy*Dz/dxw;
          ae[FIX(i,j,k)]= kapa*Dy*Dz/dxe;
          an[FIX(i,j,k)]= kapa*Dx*Dz/dyn;
          as[FIX(i,j,k)]= kapa*Dx*Dz/dys;
          af[FIX(i,j,k)]= kapa*Dx*Dy/dzf;
          ab[FIX(i,j,k)]= kapa*Dx*Dy/dzb;
        }
        
        ap0[FIX(i,j,k)]= Dx*Dy*Dz/dt;
        b[FIX(i,j,k)]= psi0[FIX(i,j,k)]*ap0[FIX(i,j,k)];
        
        if(para->prob->gravdir==GRAVZ || para->prob->gravdir==GRAVZN ) {
          b[FIX(i,j,k)] -= beta*gravz*(0.5f*(Temp[FIX(i,j,k)]
                                     +Temp[FIX(i,j,k+1)])-Temp_opt)*Dx*Dy*Dz;
          if(para->prob->plume_mod==1)
            b[FIX(i,j,k)] += s_plume[FIX(i,j,k)]/dt*Dx*Dy*Dz;//plume model 
        }
        ap[FIX(i,j,k)] = ap0[FIX(i,j,k)] + ae[FIX(i,j,k)] + aw[FIX(i,j,k)] 
                       + an[FIX(i,j,k)]  + as[FIX(i,j,k)] + af[FIX(i,j,k)]
                       + ab[FIX(i,j,k)];
      END_FOR
        
      set_bnd(para, var, var_type, psi,BINDEX);
      
      break;

    //Coefficients for temperature                           
    case TEMP: 
      kapa = para->prob->nu; 
      FOR_EACH_CELL
        if(flagp[FIX(i,j,k)]>0) continue;
        dxe=x[FIX(i+1,j,k)]-x[FIX(i,j,k)];
        dxw=x[FIX(i,j,k)]-x[FIX(i-1,j,k)];
        dyn=y[FIX(i,j+1,k)]-y[FIX(i,j,k)];
        dys=y[FIX(i,j,k)]-y[FIX(i,j-1,k)];
        dzf=z[FIX(i,j,k+1)]-z[FIX(i,j,k)];
        dzb=z[FIX(i,j,k)]-z[FIX(i,j,k-1)];
        Dx=gx[FIX(i,j,k)]-gx[FIX(i-1,j,k)];
        Dy=gy[FIX(i,j,k)]-gy[FIX(i,j-1,k)];
        Dz=gz[FIX(i,j,k)]-gz[FIX(i,j,k-1)];

        if(para->prob->tur_model == CHEN)  {
          aw[FIX(i,j,k)]= (kapa+0.5f*(vt[FIX(i,j,k)]
                                     +vt[FIX(i-1,j,k)]))/Prt*Dy*Dz/dxw;
          ae[FIX(i,j,k)]= (kapa+0.5f*(vt[FIX(i,j,k)]
                                     +vt[FIX(i+1,j,k)]))/Prt*Dy*Dz/dxe;
          an[FIX(i,j,k)]= (kapa+0.5f*(vt[FIX(i,j,k)]
                                     +vt[FIX(i,j+1,k)]))/Prt*Dx*Dz/dyn;
          as[FIX(i,j,k)]= (kapa+0.5f*(vt[FIX(i,j,k)]
                                     +vt[FIX(i,j-1,k)]))/Prt*Dx*Dz/dys;
          af[FIX(i,j,k)]= (kapa+0.5f*(vt[FIX(i,j,k)]
                                     +vt[FIX(i,j,k+1)]))/Prt*Dx*Dy/dzf;
          ab[FIX(i,j,k)]= (kapa+0.5f*(vt[FIX(i,j,k)]
                                     +vt[FIX(i,j,k-1)]))/Prt*Dx*Dy/dzb;
        }
        else {
          aw[FIX(i,j,k)]= kapa/Prt*Dy*Dz/dxw;
          ae[FIX(i,j,k)]= kapa/Prt*Dy*Dz/dxe;
          an[FIX(i,j,k)]= kapa/Prt*Dx*Dz/dyn;
          as[FIX(i,j,k)]= kapa/Prt*Dx*Dz/dys;
          af[FIX(i,j,k)]= kapa/Prt*Dx*Dy/dzf;
          ab[FIX(i,j,k)]= kapa/Prt*Dx*Dy/dzb; 
        }

        if(para->prob->plume_mod==1){

          if(BINDEX[20][FIX(i-1,j,k)]==1) aw[FIX(i,j,k)]= 0;
          if(BINDEX[20][FIX(i+1,j,k)]==1) ae[FIX(i,j,k)]= 0;
          if(BINDEX[20][FIX(i,j,k+1)]==1) af[FIX(i,j,k)]= 0;
          if(BINDEX[20][FIX(i,j,k-1)]==1) ab[FIX(i,j,k)]= 0; 
          if(BINDEX[20][FIX(i,j+1,k)]==1) an[FIX(i,j,k)]= 0;
          if(BINDEX[20][FIX(i,j-1,k)]==1) as[FIX(i,j,k)]= 0; 
        }

        ap0[FIX(i,j,k)]= Dx*Dy*Dz/dt;
        b[FIX(i,j,k)]= psi0[FIX(i,j,k)]*ap0[FIX(i,j,k)];

      END_FOR

      set_bnd(para, var, var_type, psi,BINDEX);

      FOR_EACH_CELL
        ap[FIX(i,j,k)] = ap0[FIX(i,j,k)] + ae[FIX(i,j,k)] + aw[FIX(i,j,k)] 
                       + an[FIX(i,j,k)]  + as[FIX(i,j,k)] + af[FIX(i,j,k)] 
                       + ab[FIX(i,j,k)];
      END_FOR
      break;

    //Coefficients for species    
    case DEN:  
      kapa = para->prob->nu; 
      FOR_EACH_CELL
        if(flagp[FIX(i,j,k)]>0) continue;
        dxe=x[FIX(i+1,j,k)]-x[FIX(i,j,k)];
        dxw=x[FIX(i,j,k)]-x[FIX(i-1,j,k)];
        dyn=y[FIX(i,j+1,k)]-y[FIX(i,j,k)];
        dys=y[FIX(i,j,k)]-y[FIX(i,j-1,k)];
        dzf=z[FIX(i,j,k+1)]-z[FIX(i,j,k)];
        dzb=z[FIX(i,j,k)]-z[FIX(i,j,k-1)];
        Dx=gx[FIX(i,j,k)]-gx[FIX(i-1,j,k)];
        Dy=gy[FIX(i,j,k)]-gy[FIX(i,j-1,k)];
        Dz=gz[FIX(i,j,k)]-gz[FIX(i,j,k-1)];

        if(para->prob->tur_model == CHEN)  {
          aw[FIX(i,j,k)]= (kapa+0.5f*(vt[FIX(i,j,k)]
                                     +vt[FIX(i-1,j,k)]))/Prt*Dy*Dz/dxw;
          ae[FIX(i,j,k)]= (kapa+0.5f*(vt[FIX(i,j,k)]
                                     +vt[FIX(i+1,j,k)]))/Prt*Dy*Dz/dxe;
          an[FIX(i,j,k)]= (kapa+0.5f*(vt[FIX(i,j,k)]
                                     +vt[FIX(i,j+1,k)]))/Prt*Dx*Dz/dyn;
          as[FIX(i,j,k)]= (kapa+0.5f*(vt[FIX(i,j,k)]
                                     +vt[FIX(i,j-1,k)]))/Prt*Dx*Dz/dys;
          af[FIX(i,j,k)]= (kapa+0.5f*(vt[FIX(i,j,k)]
                                     +vt[FIX(i,j,k+1)]))/Prt*Dx*Dy/dzf;
          ab[FIX(i,j,k)]= (kapa+0.5f*(vt[FIX(i,j,k)]
                                     +vt[FIX(i,j,k-1)]))/Prt*Dx*Dy/dzb;
        }
        else {
          aw[FIX(i,j,k)]= kapa/Prt*Dy*Dz/dxw;
          ae[FIX(i,j,k)]= kapa/Prt*Dy*Dz/dxe;
          an[FIX(i,j,k)]= kapa/Prt*Dx*Dz/dyn;
          as[FIX(i,j,k)]= kapa/Prt*Dx*Dz/dys;
          af[FIX(i,j,k)]= kapa/Prt*Dx*Dy/dzf;
          ab[FIX(i,j,k)]= kapa/Prt*Dx*Dy/dzb; 
        }

        ap0[FIX(i,j,k)]= Dx*Dy*Dz/dt;
        b[FIX(i,j,k)]= psi0[FIX(i,j,k)]*ap0[FIX(i,j,k)];

      END_FOR

      set_bnd(para, var, var_type, psi,BINDEX);

      FOR_EACH_CELL
        ap[FIX(i,j,k)] = ap0[FIX(i,j,k)] + ae[FIX(i,j,k)] + aw[FIX(i,j,k)] 
                       + an[FIX(i,j,k)]  + as[FIX(i,j,k)] + af[FIX(i,j,k)] 
                       + ab[FIX(i,j,k)];
      END_FOR
      break;
  }

   
}// End of coef_diff( )


