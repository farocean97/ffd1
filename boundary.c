///////////////////////////////////////////////////////////////////////////////
///
/// \file   bounary.c
///
/// \brief  Subroutines for setting boundary conditions
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
/// This file provides functions that used for the setting boundary for FFD.
/// FFD will set up bounary conditons for velocity \c set_vel_bnd(), pressure
/// \c set_bnd_pressure(),temperature \c set_bnd_temp() and species 
/// concentration \c set_bnd_density(). In addition, this fine contains 
/// functions for enforcing mass conservation at the outflow bounary and 
/// special treatment for plume cells.
///
///////////////////////////////////////////////////////////////////////////////


#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "data_structure.h"
#include "boundary.h"
#include "utility.h"




///////////////////////////////////////////////////////////////////////////////
/// \brief Entrance of setting boundary conditions
///
/// Specific boundary conditions will be selected according to the variable 
/// type.
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param var_type The type of variable
///\param psi Pointer to the variable needing the boundary conditions
///\param BINDEX Pointer to boundary index
///
///\return void No return needed
///////////////////////////////////////////////////////////////////////////////
void set_bnd(PARA_DATA *para, REAL **var, int var_type, REAL *psi, int **BINDEX){
  switch(var_type){
    case VX:
      set_bnd_vel(para, var, VX, psi, BINDEX); break;
    case VY:
      set_bnd_vel(para, var, VY, psi, BINDEX); break;
    case VZ:
      set_bnd_vel(para, var, VZ, psi, BINDEX); break;
    case TEMP:
      set_bnd_temp(para, var, TEMP, psi,BINDEX); break;
    case DEN:
      set_bnd_density(para, var,DEN,psi,BINDEX); break;
  }
} // set_bnd() 


///////////////////////////////////////////////////////////////////////////////
///\brief Set boundary conditions for velocity
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param var_type The type of variable
///\param psi Pointer to the variable needing the boundary conditions
///\param BINDEX Pointer to boundary index
///
///\return void No return needed
///////////////////////////////////////////////////////////////////////////////
void set_bnd_vel(PARA_DATA *para, REAL **var, int var_type, REAL *psi, int **BINDEX){
  int i, j, k;
  int it;
  int ii,jj,kk;
  int imax = para->geom->imax, jmax = para->geom->jmax;
  int kmax = para->geom->kmax;
  int pindexmax_x,pindexmax_y,pindexmax_z,indexmax,indexmax_us;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  REAL *aw = var[AW], *ae = var[AE], *as = var[AS], *an = var[AN];
  REAL *af = var[AF], *ab = var[AB],*ap0=var[AP0], *b=var[B],*ap=var[AP];
  REAL *gx=var[GX],*gy=var[GY],*gz=var[GZ];
  REAL *flagp = var[FLAGP],*flagu = var[FLAGU],*flagv = var[FLAGV],*flagw = var[FLAGW];


  indexmax=para->geom->index[1];
  indexmax_us=para->geom->index[2];
  pindexmax_x=para->geom->index[3];
  pindexmax_y=para->geom->index[4];
  pindexmax_z=para->geom->index[5];
  
  switch(var_type){
    //Assigning bounary condition for velocity at X direction.
    case VX:
      for(it=1;it<=indexmax;it++){
        i=BINDEX[0][it];
        j=BINDEX[1][it];
        k=BINDEX[2][it]; 
        // flagp=0 is heat source cell, then skip assigning boundary condition
        // first asssigning bounary condtion at East face of a boundary cell.
        if(flagp[FIX(i,j,k)]==0) continue;
        //flagu=0 is inlet boundary
        if(flagu[FIX(i,j,k)]==0){
          psi[FIX(i,j,k)]= para->bc->u_bc[BINDEX[13][FIX(i,j,k)]]; 
          if(fabs(psi[FIX(i,j,k)]>0.000001)) {
            b[FIX(i+1,j,k)] += para->bc->um_bc[BINDEX[13][FIX(i,j,k)]]*
                  (gy[FIX(i,j,k)]-gy[FIX(i,j-1,k)])*(gz[FIX(i,j,k)]-gz[FIX(i,j,k-1)]);
          }
        }
        //flagu>1 velocity boudnary for solid wall 
        if(flagu[FIX(i,j,k)]>1) psi[FIX(i,j,k)]= 0;

        //flag
        if(flagu[FIX(i,j,k)]==1 && i<=imax) {
          psi[FIX(i,j,k)]= psi[FIX(i+1,j,k)]; 
          ap[FIX(i+1,j,k)] -= aw[FIX(i+1,j,k)];
          aw[FIX(i+1,j,k)]=0; 
        }
        
        // Then assigning boundary condition at West face of a boundary cell
        ii=max(i-1,0);
        if(flagu[FIX(ii,j,k)]==0){
          psi[FIX(ii,j,k)]= para->bc->u_bc[BINDEX[13][FIX(ii,j,k)]]; 
          if(fabs(psi[FIX(ii,j,k)]>0.000001)) {
            b[FIX(ii-1,j,k)] += para->bc->um_bc[BINDEX[13][FIX(ii,j,k)]]*
                  (gy[FIX(i,j,k)]-gy[FIX(i,j-1,k)])*(gz[FIX(i,j,k)]-gz[FIX(i,j,k-1)]);
          }
        }
        if(flagu[FIX(ii,j,k)]>1 ) psi[FIX(ii,j,k)]= 0;
        if(flagu[FIX(ii,j,k)]==1 && ii>0){
          psi[FIX(ii,j,k)]= psi[FIX(ii-1,j,k)]; 
          ap[FIX(ii-1,j,k)] -= ae[FIX(ii-1,j,k)]; 
          ae[FIX(ii-1,j,k)]=0; 
        }
      }

      //Searching the boundary types of partition walls
      //Set velocity to zero if the normal direction of the wall is at X direction
      for(it=1;it<=pindexmax_x;it++){
        i=BINDEX[17][it];
        j=BINDEX[18][it];
        k=BINDEX[19][it];
        psi[FIX(i,j,k)]=0;
      }
      // if the normal direction of the wall is at Y or Z direction
      // set corresponding coefficient to zero
      for(it=pindexmax_x+1;it<=pindexmax_y;it++) {
        i=BINDEX[17][it];
        j=BINDEX[18][it];
        k=BINDEX[19][it];
        an[FIX(i,j,k)]=0;
        an[FIX(i-1,j,k)]=0;
        as[FIX(i,j+1,k)]=0;
        as[FIX(i-1,j+1,k)]=0;
      }
      for(it=pindexmax_y+1;it<=pindexmax_z;it++) {
        i=BINDEX[17][it];
        j=BINDEX[18][it];
        k=BINDEX[19][it];
        af[FIX(i,j,k)]=0;
        af[FIX(i-1,j,k)]=0;
        ab[FIX(i,j,k+1)]=0;
        ab[FIX(i-1,j,k+1)]=0;
      }
      break;
    //Assigning bounary condition for velocity at Y direction.   
    case VY:
      for(it=1;it<=indexmax;it++) {
        i=BINDEX[0][it];
        j=BINDEX[1][it];
        k=BINDEX[2][it];
        if(flagp[FIX(i,j,k)]==0) continue; //heat source cell 
        if(flagv[FIX(i,j,k)]==0) {
          psi[FIX(i,j,k)]= para->bc->v_bc[BINDEX[14][FIX(i,j,k)]]; 
          if(fabs(psi[FIX(i,j,k)]>0.000001)) {
            b[FIX(i,j+1,k)] += para->bc->um_bc[BINDEX[14][FIX(i,j,k)]]*
                  (gx[FIX(i,j,k)]-gx[FIX(i-1,j,k)])*(gz[FIX(i,j,k)]-gz[FIX(i,j,k-1)]);
          }
        }
        if(flagv[FIX(i,j,k)]> 1) psi[FIX(i,j,k)]= 0;
        if(flagv[FIX(i,j,k)]==1 && j<=jmax) {
          psi[FIX(i,j,k)]= psi[FIX(i,j+1,k)]; 
          ap[FIX(i,j+1,k)] -= as[FIX(i,j+1,k)]; 
          as[FIX(i,j+1,k)]=0; 
        }

        jj=max(j-1,0);
        if(flagv[FIX(i,jj,k)]==0) {
          psi[FIX(i,jj,k)]= para->bc->v_bc[BINDEX[14][FIX(i,jj,k)]];
          if(fabs(psi[FIX(i,jj,k)]>0.000001)) {
            b[FIX(i,jj-1,k)] += para->bc->um_bc[BINDEX[14][FIX(i,jj,k)]]*
                  (gx[FIX(i,j,k)]-gx[FIX(i-1,j,k)])*(gz[FIX(i,j,k)]-gz[FIX(i,j,k-1)]);
          }
        }
        if(flagv[FIX(i,jj,k)]>1) psi[FIX(i,jj,k)]= 0;
        if(flagv[FIX(i,jj,k)]==1 && jj>0){
          psi[FIX(i,jj,k)]= psi[FIX(i,jj-1,k)];
          ap[FIX(i,jj-1,k)] -= an[FIX(i,jj-1,k)];
          an[FIX(i,jj-1,k)]=0; 
        }
      }
      for(it=1;it<=pindexmax_x;it++){   
        i=BINDEX[17][it];
        j=BINDEX[18][it];
        k=BINDEX[19][it];
        ae[FIX(i,j,k)]=0;
        ae[FIX(i,j-1,k)]=0;
        aw[FIX(i+1,j,k)]=0;
        aw[FIX(i+1,j-1,k)]=0;
      }
      for(it=pindexmax_x+1;it<=pindexmax_y;it++) {
        i=BINDEX[17][it];
        j=BINDEX[18][it];
        k=BINDEX[19][it];
        psi[FIX(i,j,k)]=0;
      }
      for(it=pindexmax_y+1;it<=pindexmax_z;it++){
        i=BINDEX[17][it];
        j=BINDEX[18][it];
        k=BINDEX[19][it];
        af[FIX(i,j,k)]=0;
        af[FIX(i,j-1,k)]=0;
        ab[FIX(i,j,k+1)]=0;
        ab[FIX(i,j-1,k+1)]=0;
      }

      break;

    //Assigning bounary condition for velocity at X direction.        
    case VZ:
      for(it=1;it<=indexmax;it++){
        i=BINDEX[0][it];
        j=BINDEX[1][it];
        k=BINDEX[2][it];
        if(flagp[FIX(i,j,k)]==0) continue;
        if(flagw[FIX(i,j,k)]==0) {
          psi[FIX(i,j,k)]= para->bc->w_bc[BINDEX[15][FIX(i,j,k)]]; 
          if(fabs(psi[FIX(i,j,k)]>0.000001)) {
            b[FIX(i,j,k+1)] += para->bc->um_bc[BINDEX[15][FIX(i,j,k)]]*
                  (gx[FIX(i,j,k)]-gx[FIX(i-1,j,k)])*(gy[FIX(i,j,k)]-gy[FIX(i,j-1,k)]);
          }
        }
        if(flagw[FIX(i,j,k)]>1) psi[FIX(i,j,k)]= 0; 
        if(flagw[FIX(i,j,k)]==1 && k<=kmax){
          psi[FIX(i,j,k)]= psi[FIX(i,j,k+1)];
          ap[FIX(i,j,k+1)] -= ab[FIX(i,j,k+1)];
          ab[FIX(i,j,k+1)]=0; 
        }

        kk=max(k-1,0);
        if(flagw[FIX(i,j,kk)]==0) {
          psi[FIX(i,j,kk)]= para->bc->w_bc[BINDEX[15][FIX(i,j,kk)]]; 
          if(fabs(psi[FIX(i,j,kk)]>0.000001)) {
            b[FIX(i,j,kk-1)] += para->bc->um_bc[BINDEX[15][FIX(i,j,kk)]]*
                  (gx[FIX(i,j,k)]-gx[FIX(i-1,j,k)])*(gy[FIX(i,j,k)]-gy[FIX(i,j-1,k)]);
          }
        }
        if(flagw[FIX(i,j,kk)]>1) psi[FIX(i,j,kk)]= 0; 
        if(flagw[FIX(i,j,kk)]==1 && kk>0){
          psi[FIX(i,j,kk)]= psi[FIX(i,j,kk-1)]; 
          ap[FIX(i,j,kk-1)]-= af[FIX(i,j,kk-1)]; 
          af[FIX(i,j,kk-1)]=0; 
        }
      }
      for(it=1;it<=pindexmax_x;it++) {
        i=BINDEX[17][it];
        j=BINDEX[18][it];
        k=BINDEX[19][it];
        ae[FIX(i,j,k)]=0;
        ae[FIX(i,j,k-1)]=0;
        aw[FIX(i+1,j,k)]=0;
        aw[FIX(i+1,j,k-1)]=0;
      }
      for(it=pindexmax_x+1;it<=pindexmax_y;it++) {
        i=BINDEX[17][it];
        j=BINDEX[18][it];
        k=BINDEX[19][it];
        an[FIX(i,j,k)]=0;
        an[FIX(i,j,k-1)]=0;
        as[FIX(i,j+1,k)]=0;
        as[FIX(i,j+1,k-1)]=0;
      }
      for(it=pindexmax_y+1;it<=pindexmax_z;it++) {
        i=BINDEX[17][it];
        j=BINDEX[18][it];
        k=BINDEX[19][it];
        psi[FIX(i,j,k)]=0;
      }
      break;
   }
}// End of set_bnd_vel( )




///////////////////////////////////////////////////////////////////////////////
///\brief Set boundary conditions for temperature
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param var_type The type of variable
///\param psi Pointer to the variable needing the boundary conditions
///\param BINDEX Pointer to boundary index
///
///\return void No return needed
///////////////////////////////////////////////////////////////////////////////

void set_bnd_temp(PARA_DATA *para, REAL **var, int var_type, REAL *psi,int **BINDEX){
  int i, j, k;
  int it,zone_num;
  int pindexmax_x,pindexmax_y,pindexmax_z,indexmax,indexmax_us;
  int imax = para->geom->imax, jmax = para->geom->jmax;
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  REAL *aw = var[AW], *ae = var[AE], *as = var[AS], *an = var[AN];
  REAL *af = var[AF], *ab = var[AB],*ap=var[AP],*b=var[B],*q=var[FLUX];
  REAL *gx = var[GX], *gy = var[GY], *gz = var[GZ];
  REAL *fx = var[FX], *fy = var[FY], *fz = var[FZ];
  REAL *T=var[TEMP];
  REAL coeff_h=para->prob->coeff_h;
  REAL rho=para->prob->rho;
  REAL cp=para->prob->spec;
  REAL T_amb=292.5f;
  REAL qflux;
  REAL qflow_a=para->bc->qflow_a;
  int priority=0;//the priority of storing temperature value at the cell.
  REAL *flagp = var[FLAGP],*flagu = var[FLAGU];
  REAL *flagv = var[FLAGV],*flagw = var[FLAGW];

  indexmax=para->geom->index[1];
  indexmax_us=para->geom->index[2];
  pindexmax_x=para->geom->index[3];
  pindexmax_y=para->geom->index[4];
  pindexmax_z=para->geom->index[5];

  for(it=1;it<=indexmax;it++){
    priority=0;
    i=BINDEX[0][it];
    j=BINDEX[1][it];
    k=BINDEX[2][it];
    
    //X direction EAST FACE
    if(i<imax && (flagp[FIX(i+1,j,k)]<=0)) {
      switch((int)flagu[FIX(i,j,k)]){
        case 0: //inlet_surface
          zone_num=BINDEX[13][FIX(i,j,k)];
          priority=4;
          psi[FIX(i,j,k)]=para->bc->t_bc[zone_num];
          break;
        case 1: //outlet_surface
          aw[FIX(i+1,j,k)]=0;
          if(priority<3){
            psi[FIX(i,j,k)]=psi[FIX(i+1,j,k)];
            priority=3;
          }
          break;
        case 2: //wall_surface
          zone_num=BINDEX[13][FIX(i,j,k)];
          if(para->bc->fltmp[zone_num]==0) {
            aw[FIX(i+1,j,k)]=0;
            qflux=para->bc->t_bc[zone_num];
            b[FIX(i+1,j,k)] += 1/(rho*cp)*qflux*(gy[FIX(i,j,k)]-gy[FIX(i,j-1,k)])
                                               *(gz[FIX(i,j,k)]-gz[FIX(i,j,k-1)]);
            if(priority<2){
              psi[FIX(i,j,k)]= qflux/4.0f+psi[FIX(i+1,j,k)];
              priority=2;
            }  
          }
          else {
            aw[FIX(i+1,j,k)]= coeff_h*(gy[FIX(i,j,k)]-gy[FIX(i,j-1,k)])
                                     *(gz[FIX(i,j,k)]-gz[FIX(i,j,k-1)]);
            if(priority<2){
              psi[FIX(i,j,k)]= para->bc->t_bc[zone_num];
              priority=2;
            }
          }
          break;
        case 3://block_wall_surface
          zone_num=BINDEX[16][FIX(i,j,k)];
          if(para->bc->fltmp[zone_num]==0) {
            aw[FIX(i+1,j,k)]=0;
            qflux=para->bc->t_bc[zone_num];
            b[FIX(i+1,j,k)]+=1/(rho*cp)*qflux*(gy[FIX(i,j,k)]-gy[FIX(i,j-1,k)])
                                             *(gz[FIX(i,j,k)]-gz[FIX(i,j,k-1)]);
            if(priority<=1){
              psi[FIX(i,j,k)]= qflux/4.0f+psi[FIX(i+1,j,k)]; 
              priority=1;
            }  
          }
          else {
            aw[FIX(i+1,j,k)]= coeff_h*(gy[FIX(i,j,k)]-gy[FIX(i,j-1,k)])
                                     *(gz[FIX(i,j,k)]-gz[FIX(i,j,k-1)]);
            if(priority<=1){
              psi[FIX(i,j,k)]= para->bc->t_bc[zone_num]; 
              priority=1;
            } 
          }
          break;
      }
    }
    
    //X direction WEST FACE
    if(i>0 && (flagp[FIX(i-1,j,k)]<=0)) {
      switch((int)flagu[FIX(i-1,j,k)]){
        case 0:
          zone_num=BINDEX[13][FIX(i-1,j,k)];
          priority=4;
          psi[FIX(i,j,k)]=para->bc->t_bc[zone_num];
          break;
        case 1:
          ae[FIX(i-1,j,k)]=0;
          if(priority<3){
            psi[FIX(i,j,k)]=psi[FIX(i-1,j,k)];
            priority=3;
          }
          break;
        case 2:
          zone_num=BINDEX[13][FIX(i-1,j,k)];
          if(para->bc->fltmp[zone_num]==0) {
            ae[FIX(i-1,j,k)]=0;
            qflux=para->bc->t_bc[zone_num];
            b[FIX(i-1,j,k)] += 1/(rho*cp)*qflux*(gy[FIX(i,j,k)]-gy[FIX(i,j-1,k)])
                                             *(gz[FIX(i,j,k)]-gz[FIX(i,j,k-1)]);
            if(priority<2){
              psi[FIX(i,j,k)]= qflux/4.0f+psi[FIX(i-1,j,k)];
              priority=2;
            } 
          }
          else {
            ae[FIX(i-1,j,k)]= coeff_h*(gy[FIX(i,j,k)]-gy[FIX(i,j-1,k)])
                                     *(gz[FIX(i,j,k)]-gz[FIX(i,j,k-1)]);
            if(priority<2){
              psi[FIX(i,j,k)]= para->bc->t_bc[zone_num];
              priority=2;
            } 
          }
          break;
        case 3:
          zone_num=BINDEX[16][FIX(i,j,k)];
          if(para->bc->fltmp[zone_num]==0){
            ae[FIX(i-1,j,k)]=0;
            qflux=para->bc->t_bc[zone_num];
            b[FIX(i-1,j,k)] += 1/(rho*cp)*qflux*(gy[FIX(i,j,k)]-gy[FIX(i,j-1,k)])
                                               *(gz[FIX(i,j,k)]-gz[FIX(i,j,k-1)]);
            if(priority<=1) {
              psi[FIX(i,j,k)]= qflux/4.0f+psi[FIX(i-1,j,k)]; 
              priority=1;
            } 
          }
          else{
            ae[FIX(i-1,j,k)]= coeff_h*(gy[FIX(i,j,k)]-gy[FIX(i,j-1,k)])
                                     *(gz[FIX(i,j,k)]-gz[FIX(i,j,k-1)]);
            if(priority<=1){
              psi[FIX(i,j,k)]= para->bc->t_bc[zone_num];
              priority=1;
            }
          }
          break;
      }
    }
    
    //Y direction NORTH FACE
    if(j<jmax && (flagp[FIX(i,j+1,k)]<=0)) {
      switch( (int)flagv[FIX(i,j,k)]){
        case 0:
          zone_num=BINDEX[14][FIX(i,j,k)];
          priority=4;
          psi[FIX(i,j,k)]=para->bc->t_bc[zone_num];
          break;
        case 1:
          as[FIX(i,j+1,k)]=0;
          if(priority<3){
            psi[FIX(i,j,k)]=psi[FIX(i,j+1,k)];
            priority=3;
          }
          break;
        case 2:
          zone_num=BINDEX[14][FIX(i,j,k)];
          if(para->bc->fltmp[zone_num]==0){
            as[FIX(i,j+1,k)]=0;
            qflux=para->bc->t_bc[zone_num];
            b[FIX(i,j+1,k)]+=1/(rho*cp)*qflux*(gx[FIX(i,j,k)]-gx[FIX(i-1,j,k)])
                                             *(gz[FIX(i,j,k)]-gz[FIX(i,j,k-1)]);
            if(priority<2){
              psi[FIX(i,j,k)]= qflux/4.0f+psi[FIX(i,j+1,k)]; 
              priority=2;
            }
          }
          else{
            as[FIX(i,j+1,k)]= coeff_h*(gx[FIX(i,j,k)]-gx[FIX(i-1,j,k)])
                                     *(gz[FIX(i,j,k)]-gz[FIX(i,j,k-1)]);
            if(priority<2){
              psi[FIX(i,j,k)]= para->bc->t_bc[zone_num]; 
              priority=2;
            }
          }
          break;
        case 3:
          zone_num=BINDEX[16][FIX(i,j,k)];
          if(para->bc->fltmp[zone_num]==0) {
            as[FIX(i,j+1,k)]=0;
            qflux=para->bc->t_bc[zone_num];
            b[FIX(i,j+1,k)]+=1/(rho*cp)*qflux*(gx[FIX(i,j,k)]-gx[FIX(i-1,j,k)])
                                             *(gz[FIX(i,j,k)]-gz[FIX(i,j,k-1)]);
            if(priority<=1){
              psi[FIX(i,j,k)]= qflux/4.0f+psi[FIX(i,j+1,k)]; 
              priority=1;
            }
          }
          else {
            as[FIX(i,j+1,k)]= coeff_h*(gx[FIX(i,j,k)]-gx[FIX(i-1,j,k)])
                                     *(gz[FIX(i,j,k)]-gz[FIX(i,j,k-1)]);
            if(priority<=1){
              psi[FIX(i,j,k)]= para->bc->t_bc[zone_num];
              priority=1;
            } 
          }
          break;
      }
    }
    
    
    //Y direction SOUTH FACE
    if(j>0 && (flagp[FIX(i,j-1,k)]<=0)){
      switch( (int)flagv[FIX(i,j-1,k)]){
        case 0:
          zone_num=BINDEX[14][FIX(i,j-1,k)];
          priority=4;
          psi[FIX(i,j,k)]=para->bc->t_bc[zone_num];
          break;
        case 1:
          an[FIX(i,j-1,k)]=0;
          if(priority<3){
            psi[FIX(i,j,k)]=psi[FIX(i,j-1,k)];
            priority=3;
          }
          break;
        case 2:
          zone_num=BINDEX[14][FIX(i,j-1,k)];
          if(para->bc->fltmp[zone_num]==0){
            an[FIX(i,j-1,k)]=0;
            qflux=para->bc->t_bc[zone_num];
            b[FIX(i,j-1,k)]+=1/(rho*cp)*qflux*(gx[FIX(i,j,k)]-gx[FIX(i-1,j,k)])
                                             *(gz[FIX(i,j,k)]-gz[FIX(i,j,k-1)]);
            if(priority<2){
              psi[FIX(i,j,k)]= qflux/4.0f+psi[FIX(i,j-1,k)];
              priority=2;
            }
          }
          else{
            an[FIX(i,j-1,k)]= coeff_h*(gx[FIX(i,j,k)]-gx[FIX(i-1,j,k)])
                                     *(gz[FIX(i,j,k)]-gz[FIX(i,j,k-1)]);
            if(priority<2) {
              psi[FIX(i,j,k)]= para->bc->t_bc[zone_num];
              priority=2;
            } 
          }
          break;
        case 3:
          zone_num=BINDEX[16][FIX(i,j,k)];
          if(para->bc->fltmp[zone_num]==0){
            an[FIX(i,j-1,k)]=0;
            qflux=para->bc->t_bc[zone_num];
            b[FIX(i,j-1,k)] += 1/(rho*cp)*qflux*(gx[FIX(i,j,k)]-gx[FIX(i-1,j,k)])
                                               *(gz[FIX(i,j,k)]-gz[FIX(i,j,k-1)]);
            if(priority<=1){
              psi[FIX(i,j,k)]= qflux/4.0f+psi[FIX(i,j-1,k)];
              priority=1;
            }
          }
          else{
            an[FIX(i,j-1,k)]= coeff_h*(gx[FIX(i,j,k)]-gx[FIX(i-1,j,k)])
                                     *(gz[FIX(i,j,k)]-gz[FIX(i,j,k-1)]);
            if(priority<=1) {
              psi[FIX(i,j,k)]= para->bc->t_bc[zone_num];
              priority=1;
            }
          }
          break;
      }
    }


    //Z direction FRONT FACE
    if(k<kmax && (flagp[FIX(i,j,k+1)]<=0)){
      switch( (int)flagw[FIX(i,j,k)]){
        case 0:
          zone_num=BINDEX[15][FIX(i,j,k)];
          priority=4;
          psi[FIX(i,j,k)]=para->bc->t_bc[zone_num];
          break;
        case 1:
          ab[FIX(i,j,k+1)]=0;
          if(priority<3){
            psi[FIX(i,j,k)]=psi[FIX(i,j,k+1)];
            priority=3;
          }
          break;
        case 2:
          zone_num=BINDEX[15][FIX(i,j,k)];
          if(para->bc->fltmp[zone_num]==0){
            ab[FIX(i,j,k+1)]=0;
            qflux=para->bc->t_bc[zone_num];
            b[FIX(i,j,k+1)]+=1/(rho*cp)*qflux*(gy[FIX(i,j,k)]-gy[FIX(i,j-1,k)])
                                             *(gx[FIX(i,j,k)]-gx[FIX(i-1,j,k)]);
            if(priority<2){
              psi[FIX(i,j,k)]= qflux/4.0f+psi[FIX(i,j,k+1)];
              priority=2;
            }
          }
          else{
            ab[FIX(i,j,k+1)]= coeff_h*(gy[FIX(i,j,k)]-gy[FIX(i,j-1,k)])
                                     *(gx[FIX(i,j,k)]-gx[FIX(i-1,j,k)]);
            if(priority<2){
              psi[FIX(i,j,k)]= para->bc->t_bc[zone_num];
              priority=2;
            } 
          }
          break;
        case 3:
          zone_num=BINDEX[16][FIX(i,j,k)];
          if(para->bc->fltmp[zone_num]==0){
            ab[FIX(i,j,k+1)]=0;
            qflux=para->bc->t_bc[zone_num];
            b[FIX(i,j,k+1)] += 1/(rho*cp)*qflux*(gy[FIX(i,j,k)]-gy[FIX(i,j-1,k)])
                                               *(gx[FIX(i,j,k)]-gx[FIX(i-1,j,k)]);
            if(priority<=1){
              psi[FIX(i,j,k)]= qflux/4.0f+psi[FIX(i,j,k+1)];
              priority=1;
            } 
          }
          else {
            ab[FIX(i,j,k+1)]= coeff_h*(gy[FIX(i,j,k)]-gy[FIX(i,j-1,k)])
                                     *(gx[FIX(i,j,k)]-gx[FIX(i-1,j,k)]);
            if(priority<=1) {
              psi[FIX(i,j,k)]= para->bc->t_bc[zone_num];
              priority=1; 
            }
          }
          break;
      }
    }
    
    
    //Z direction BACK FACE
    if(k>0 && (flagp[FIX(i,j,k-1)]<=0)){
      switch( (int)flagw[FIX(i,j,k-1)]){
        case 0:
          zone_num=BINDEX[15][FIX(i,j,k-1)];
          priority=4;
          psi[FIX(i,j,k)]=para->bc->t_bc[zone_num];
          break;
        case 1:
          af[FIX(i,j,k-1)]=0;
          if(priority<3){
            psi[FIX(i,j,k)]=psi[FIX(i,j,k-1)];
            priority=3;
          }
          break;
        case 2:
          zone_num=BINDEX[15][FIX(i,j,k-1)];
          if(para->bc->fltmp[zone_num]==0){
            af[FIX(i,j,k-1)]=0;
            qflux=para->bc->t_bc[zone_num];
            b[FIX(i,j,k-1)] += 1/(rho*cp)*qflux*(gy[FIX(i,j,k)]-gy[FIX(i,j-1,k)])
                                               *(gx[FIX(i,j,k)]-gx[FIX(i-1,j,k)]);
            if(priority<2){
              psi[FIX(i,j,k)]= qflux/4.0f+psi[FIX(i,j,k-1)]; 
              priority=2;
            }
          }
          else {
            af[FIX(i,j,k-1)]= coeff_h*(gy[FIX(i,j,k)]-gy[FIX(i,j-1,k)])
                                     *(gx[FIX(i,j,k)]-gx[FIX(i-1,j,k)]);
            if(priority<2){
              psi[FIX(i,j,k)]= para->bc->t_bc[zone_num];
              priority=2;
            }
          }
          break;
        case 3:
          zone_num=BINDEX[16][FIX(i,j,k)];
          if(para->bc->fltmp[zone_num]==0){
            af[FIX(i,j,k-1)]=0;
            qflux=para->bc->t_bc[zone_num];
            b[FIX(i,j,k-1)] += 1/(rho*cp)*qflux*(gy[FIX(i,j,k)]-gy[FIX(i,j-1,k)])
                                               *(gx[FIX(i,j,k)]-gx[FIX(i-1,j,k)]);
            if(priority<=1) {
              psi[FIX(i,j,k)]= qflux/4.0f+psi[FIX(i,j,k-1)];
              priority=1;
            }
          }
          else{
            af[FIX(i,j,k-1)]= coeff_h*(gy[FIX(i,j,k)]-gy[FIX(i,j-1,k)])
                                     *(gx[FIX(i,j,k)]-gx[FIX(i-1,j,k)]);
            if(priority<=1) {
              psi[FIX(i,j,k)]= para->bc->t_bc[zone_num];
              priority=1;
            }
          }
          break;
      }
    }
  }
  
  
  //Assigning temperature boundary conditions for partition wall
  //Assumption: these walls are adiabatic
  for(it=1;it<=pindexmax_x;it++) {
    i=BINDEX[17][it];
    j=BINDEX[18][it];
    k=BINDEX[19][it];
    ae[FIX(i,j,k)]=0;
    aw[FIX(i+1,j,k)]=0;
  }
  for(it=pindexmax_x+1;it<=pindexmax_y;it++) {
    i=BINDEX[17][it];
    j=BINDEX[18][it];
    k=BINDEX[19][it];
    an[FIX(i,j,k)]=0;
    as[FIX(i,j+1,k)]=0;
  }
  for(it=pindexmax_y+1;it<=pindexmax_z;it++) {
    i=BINDEX[17][it];
    j=BINDEX[18][it];
    k=BINDEX[19][it];
    af[FIX(i,j,k)]=0;
    ab[FIX(i,j,k+1)]=0;
  }
  
  //Assigning boundary condtions for heat source cell
  for(it=indexmax+1;it<=indexmax_us;it++) {
    i=BINDEX[0][it];
    j=BINDEX[1][it];
    k=BINDEX[2][it];
    zone_num=BINDEX[16][FIX(i,j,k)];

    b[FIX(i,j,k)] += para->bc->t_bc[zone_num]/(rho*cp);

  }
 
} //End of set_bnd_temp()


///////////////////////////////////////////////////////////////////////////////
///\brief Set boundary conditions for species
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param var_type The type of variable
///\param p Pointer to the variable needing the boundary conditions
///\param BINDEX Pointer to boundary index
///
///\return void No return needed
///////////////////////////////////////////////////////////////////////////////
void set_bnd_density(PARA_DATA *para, REAL **var, int var_type, REAL *p,int **BINDEX) {
  int i, j, k;
  int imax = para->geom->imax, jmax = para->geom->jmax;
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  REAL *aw = var[AW], *ae = var[AE], *as = var[AS], *an = var[AN];
  REAL *af = var[AF], *ab = var[AB];
  REAL *flagp = var[FLAGP];

  FOR_ALL_CELL
    if(flagp[FIX(i,j,k)]>0) continue;
    if(flagp[FIX(i,j,k)]==0) {
      p[FIX(i,j,k)]=0;
      p[FIX(0,11,11)]=100;
    }
    else{
      if(flagp[FIX(i-1,j,k)]>0){
        p[FIX(i-1,j,k)]=p[FIX(i,j,k)];
        aw[FIX(i,j,k)]=0;
      }
      if(flagp[FIX(i+1,j,k)]>0){
        p[FIX(i+1,j,k)]=p[FIX(i,j,k)]; 
        ae[FIX(i,j,k)]=0;
      }
      if(flagp[FIX(i,j-1,k)]>0)  {
        p[FIX(i,j-1,k)]=p[FIX(i,j,k)]; 
        as[FIX(i,j,k)]=0;
      }
      if(flagp[FIX(i,j+1,k)]>0)  {
        p[FIX(i,j+1,k)]=p[FIX(i,j,k)];
        an[FIX(i,j,k)]=0;
      }
      if(flagp[FIX(i,j,k-1)]>0)  {
        p[FIX(i,j,k-1)]=p[FIX(i,j,k)];  
        ab[FIX(i,j,k)]=0;
      }
      if(flagp[FIX(i,j,k+1)]>0)  {
        p[FIX(i,j,k+1)]=p[FIX(i,j,k)]; 
        af[FIX(i,j,k)]=0;
      }
    }
    END_FOR

} // End of set_bnd_pressure( )


///////////////////////////////////////////////////////////////////////////////
///\brief Set boundary conditions for pressure equation
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param p Pointer to the variable needing the boundary conditions
///\param BINDEX Pointer to boundary index
///
///\return void No return needed
///////////////////////////////////////////////////////////////////////////////
void set_bnd_pressure(PARA_DATA *para, REAL **var, REAL *p, int **BINDEX){
  int i, j, k,it;
  int imax = para->geom->imax, jmax = para->geom->jmax;
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  REAL *aw = var[AW], *ae = var[AE], *as = var[AS], *an = var[AN];
  REAL *af = var[AF], *ab = var[AB];
  int pindexmax_x,pindexmax_y,pindexmax_z;
  REAL *flagp = var[FLAGP];

  pindexmax_x=para->geom->index[3];
  pindexmax_y=para->geom->index[4];
  pindexmax_z=para->geom->index[5];

  FOR_EACH_CELL
    if(flagp[FIX(i,j,k)]==1) continue;
    if(flagp[FIX(i-1,j,k)]==1)  {
      p[FIX(i-1,j,k)]=p[FIX(i,j,k)]; 
      aw[FIX(i,j,k)]=0;
    }
    if(flagp[FIX(i+1,j,k)]==1)  {
      p[FIX(i+1,j,k)]=p[FIX(i,j,k)];    
      ae[FIX(i,j,k)]=0;
    }
    if(flagp[FIX(i,j-1,k)]==1)  {
      p[FIX(i,j-1,k)]=p[FIX(i,j,k)];   
      as[FIX(i,j,k)]=0;
    }
    if(flagp[FIX(i,j+1,k)]==1)  {
      p[FIX(i,j+1,k)]=p[FIX(i,j,k)];   
      an[FIX(i,j,k)]=0;
    }
    if(flagp[FIX(i,j,k-1)]==1)  {
      p[FIX(i,j,k-1)]=p[FIX(i,j,k)];    
      ab[FIX(i,j,k)]=0;
    }
    if(flagp[FIX(i,j,k+1)]==1)  {
      p[FIX(i,j,k+1)]=p[FIX(i,j,k)];   
      af[FIX(i,j,k)]=0;
    }
  END_FOR
    
  //Set the boundary condtions for partition walls  
  for(it=1;it<=pindexmax_x;it++) {
    i=BINDEX[17][it];
    j=BINDEX[18][it];
    k=BINDEX[19][it];
    ae[FIX(i,j,k)]=0;
    aw[FIX(i+1,j,k)]=0;
  }

  for(it=pindexmax_x+1;it<=pindexmax_y;it++) {
    i=BINDEX[17][it];
    j=BINDEX[18][it];
    k=BINDEX[19][it];
    an[FIX(i,j,k)]=0;
    as[FIX(i,j+1,k)]=0;
  }

  for(it=pindexmax_y+1;it<=pindexmax_z;it++) {
    i=BINDEX[17][it];
    j=BINDEX[18][it];
    k=BINDEX[19][it];
    af[FIX(i,j,k)]=0;
    ab[FIX(i,j,k+1)]=0;
  }

} // End of set_bnd_pressure( )



///////////////////////////////////////////////////////////////////////////////
///\brief correct the velocity at the outflow boundary
///
/// The exterpolation is applied to calcuate the velocity value at the outflow 
/// bounary ,but it can not enforce the mass conservation.
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param var_type The type of variable
///\param psi Pointer to the variable needing the boundary conditions
///\param BINDEX Pointer to boundary index
///
///\return void No return needed
///////////////////////////////////////////////////////////////////////////////
void mass_conservation(PARA_DATA *para, REAL **var,int **BINDEX) {
  int i, j, k;
  int SI,EI,SJ,EJ,SK,EK;
  int imax = para->geom->imax, jmax = para->geom->jmax;
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  REAL *gx = var[GX], *gy = var[GY], *gz = var[GZ];
  REAL *u = var[VXM], *v = var[VYM], *w = var[VZM];
  REAL mass_in=0, mass_out = 0,mass_out2=0;
  REAL area_out=0;
  REAL *flagp = var[FLAGP];
  int zone_num,zone_inlet,zone_outlet;


  zone_inlet=para->geom->zone_inlet;
  zone_outlet=para->geom->zone_outlet;
  
  
  // compute the total inflow rate
  
  for(zone_num=1;zone_num<=zone_inlet;zone_num++) {
    SI=BINDEX[7][zone_num];
    EI=BINDEX[8][zone_num];
    SJ=BINDEX[9][zone_num];
    EJ=BINDEX[10][zone_num];
    SK=BINDEX[11][zone_num];
    EK=BINDEX[12][zone_num];
    
    if(SI==EI) {
      for(j=SJ ;j<EJ ;j++)
        for(k=SK ;k<EK ;k++) {
          mass_in += flagp[FIX(SI,j+1,k+1)]*u[FIX(SI,j+1,k+1)]
                     *(gy[FIX(SI,j+1,k+1)]-gy[FIX(SI,j,k+1)])* (gz[FIX(SI,j+1,k+1)]-gz[FIX(SI,j+1,k)]);
        }
    }
    if(SJ==EJ)  {
      for(i=SI ;i<EI ;i++)
        for(k=SK ;k<EK ;k++) {
          mass_in += flagp[FIX(i+1,SJ,k+1)]*v[FIX(i+1,SJ,k+1)]
                     *(gx[FIX(i+1,SJ,k+1)]-gx[FIX(i,SJ,k+1)])* (gz[FIX(i+1,SJ,k+1)]-gz[FIX(i+1,SJ,k)]);
        }
    }
    if(SK==EK)  {
      for(i=SI ;i<EI ;i++)
        for(j=SJ ;j<EJ ;j++) {
          mass_in += flagp[FIX(i+1,j+1,SK)]*w[FIX(i+1,j+1,SK)]
                     *(gx[FIX(i+1,j+1,SK)]-gx[FIX(i,j+1,SK)])* (gy[FIX(i+1,j+1,SK)]-gy[FIX(i+1,j,SK)]);
        }
    }
  }
  
  
  
  
  //compute the total outflow rate
  
  for(zone_num=zone_inlet+1;zone_num<=zone_outlet;zone_num++) {
    SI=BINDEX[7][zone_num];
    EI=BINDEX[8][zone_num];
    SJ=BINDEX[9][zone_num];
    EJ=BINDEX[10][zone_num];
    SK=BINDEX[11][zone_num];
    EK=BINDEX[12][zone_num];
    
    
    if(SI==EI) {
      for(j=SJ ;j<EJ ;j++)
        for(k=SK ;k<EK ;k++) { 
          mass_out += flagp[FIX(SI,j+1,k+1)]*u[FIX(SI,j+1,k+1)]
                     *(gy[FIX(SI,j+1,k+1)]-gy[FIX(SI,j,k+1)])* (gz[FIX(SI,j+1,k+1)]-gz[FIX(SI,j+1,k)]);
          area_out +=(gy[FIX(SI,j+1,k+1)]-gy[FIX(SI,j,k+1)])* (gz[FIX(SI,j+1,k+1)]-gz[FIX(SI,j+1,k)]); 
        }
    }
    
    if(SJ==EJ) {
      for(i=SI ;i<EI ;i++)
        for(k=SK ;k<EK ;k++) {
          mass_out += flagp[FIX(i+1,SJ,k+1)]*v[FIX(i+1,SJ,k+1)]
                     *(gx[FIX(i+1,SJ,k+1)]-gx[FIX(i,SJ,k+1)])* (gz[FIX(i+1,SJ,k+1)]-gz[FIX(i+1,SJ,k)]);
          area_out +=(gx[FIX(i+1,SJ,k+1)]-gx[FIX(i,SJ,k+1)])* (gz[FIX(i+1,SJ,k+1)]-gz[FIX(i+1,SJ,k)]);
        }
    }
    
    if(SK==EK)  { 
      for(i=SI ;i<EI ;i++)
        for(j=SJ ;j<EJ ;j++) {
           mass_out += flagp[FIX(i+1,j+1,SK)]*w[FIX(i+1,j+1,SK)]
                     *(gx[FIX(i+1,j+1,SK)]-gx[FIX(i,j+1,SK)])* (gy[FIX(i+1,j+1,SK)]-gy[FIX(i+1,j,SK)]);
           area_out +=(gx[FIX(i+1,j+1,SK)]-gx[FIX(i,j+1,SK)])* (gy[FIX(i+1,j+1,SK)]-gy[FIX(i+1,j,SK)]);
        }
    }
  }
  

  
  // correct the velocity at the outlet
  for(zone_num=zone_inlet+1;zone_num<=zone_outlet;zone_num++) {
    SI=BINDEX[7][zone_num];
    EI=BINDEX[8][zone_num];
    SJ=BINDEX[9][zone_num];
    EJ=BINDEX[10][zone_num];
    SK=BINDEX[11][zone_num];
    EK=BINDEX[12][zone_num];
    
    
    if(SI==EI) {
      for(j=SJ ;j<EJ ;j++)
        for(k=SK ;k<EK ;k++)  {
          //u[FIX(SI,j+1,k+1)] -= flagp[FIX(SI,j+1,k+1)]*(mass_in+mass_out)/area_out; 
          u[FIX(SI,j+1,k+1)]=-flagp[FIX(SI,j+1,k+1)]*mass_in/area_out;
        }
    }
    
    if(SJ==EJ) {
      for(i=SI ;i<EI ;i++)
        for(k=SK ;k<EK ;k++)  {
          //v[FIX(i+1,SJ,k+1)] -= flagp[FIX(i+1,SJ,k+1)]*(mass_in+mass_out)/area_out;
          v[FIX(i+1,SJ,k+1)]=-flagp[FIX(i+1,SJ,k+1)]*mass_in/area_out;
        }
    }
    
    if(SK==EK) {
      for(i=SI ;i<EI ;i++)
        for(j=SJ ;j<EJ ;j++) {
          //w[FIX(i+1,j+1,SK)] -= flagp[FIX(i+1,j+1,SK)]*(mass_in+mass_out)/area_out;
          w[FIX(i+1,j+1,SK)]= -flagp[FIX(i+1,j+1,SK)]*mass_in/area_out;
        }
    }
  }

}// End of mass_conservation()



///////////////////////////////////////////////////////////////////////////////
///\brief Compute the velocity using plume model
///
/// For the coarse-grid FFD,if the heat source cell is too large,the predicted
/// velocity at the heat source cell might not be accurate, the plume model can
/// be applied to predict that veolcity, and to correct the FFD prediction.
///
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param BINDEX Pointer to boundary index
///
///\return void No return needed
///////////////////////////////////////////////////////////////////////////////

void plume(PARA_DATA *para, REAL **var, int **BINDEX){
  int i, j, k;
  int it;
  int imax = para->geom->imax, jmax = para->geom->jmax;
  int kmax = para->geom->kmax;
  int indexmax,indexmax_us;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  int zone_num;
  REAL *x = var[X],  *y = var[Y],  *z = var[Z]; 
  REAL *gx = var[GX],  *gy = var[GY],  *gz = var[GZ]; 
  REAL *u = var[VX], *v = var[VY], *w = var[VZ];
  REAL *T=var[TEMP];
  REAL *flagp = var[FLAGP],*flagv=var[FLAGV];
  REAL *vt=var[VT];
  REAL Prt=para->prob->Prt;
  REAL rho=para->prob->rho;
  REAL cp=para->prob->spec;
  REAL beta=para->prob->beta;
  REAL b;
  REAL a1,Dx,Dy,Dz;
  REAL vc,v_plume;
  REAL zv=para->geom->zv;
  REAL z_s;
  REAL Kv=0.122;
  REAL Kt=0.45; 
  REAL Kq=0.0045;
  REAL eta=0.09;
  REAL Pc;
  REAL *s_plume=var[ER];
  
  indexmax=para->geom->index[1];
  indexmax_us=para->geom->index[2];


  for(it=indexmax+1;it<=indexmax_us;it++) {
    i=BINDEX[0][it];
    j=BINDEX[1][it];
    k=BINDEX[2][it];
    zone_num=BINDEX[16][FIX(i,j,k)]; 
    Pc=para->bc->t_bc[zone_num];

    para->geom->jplume=j;
    para->geom->iplume=i;
    para->geom->kplume=k;

    switch(para->prob->gravdir) {
      case GRAVX:
      case GRAVXN:
        z_s= x[FIX(i,j,k)]+zv;
        Dx = gx[FIX(i,j,k)]-gx[FIX(i-1,j,k)];
        Dy = gy[FIX(i,j ,k)]-gy[FIX(i,j-1,k)];
        Dz = gz[FIX(i,j ,k)]-gz[FIX(i,j,k-1)];
        //plume model coefficient
        b=6/5*eta*(gx[FIX(i,j,k)]-z_s);
        a1=3.14f*b*b*(1-exp(-Dy*Dz/(b*b)))/(Dy*Dz);
        vc=Kv*pow(Pc,0.33333)*pow((gx[FIX(i,j,k)]-z_s),-0.33333);
        v_plume=a1*vc;//0.63
        //calculating the source term in momentum equation
        s_plume[FIX(i,j,k)]=v_plume-u[FIX(i,j,k)];
        break;

      case GRAVY:
      case GRAVYN:
        z_s= y[FIX(i,j,k)]+zv;
        Dx = gx[FIX(i,j,k)]-gx[FIX(i-1,j,k)];
        Dy = gy[FIX(i,j ,k)]-gy[FIX(i,j-1,k)];
        Dz = gz[FIX(i,j ,k)]-gz[FIX(i,j,k-1)];
        b=6/5*eta*(gy[FIX(i,j,k)]-z_s);
        a1=3.14f*b*b*(1-exp(-Dx*Dz/(b*b)))/(Dx*Dz);
        vc=Kv*pow(Pc,0.33333)*pow((gy[FIX(i,j,k)]-z_s),-0.33333);
        v_plume=a1*vc;//0.63
        s_plume[FIX(i,j,k)]=v_plume-v[FIX(i,j,k)];
        break;

      case GRAVZ:
      case GRAVZN:

        z_s= z[FIX(i,j,k)]+zv;
        Dx = gx[FIX(i,j,k)]-gx[FIX(i-1,j,k)];
        Dy = gy[FIX(i,j ,k)]-gy[FIX(i,j-1,k)];
        Dz = gz[FIX(i,j ,k)]-gz[FIX(i,j,k-1)];
        b=6/5*eta*(gz[FIX(i,j,k)]-z_s);
        a1=3.14f*b*b*(1-exp(-Dx*Dy/(b*b)))/(Dx*Dy);
        vc=Kv*pow(Pc,0.33333)*pow((gz[FIX(i,j,k)]-z_s),-0.33333);
        v_plume=a1*vc;//0.63
        s_plume[FIX(i,j,k)]=v_plume-w[FIX(i,j,k)];
        break;
    
    }
    
  }
  
}


///////////////////////////////////////////////////////////////////////////////
///\brief Compute the heat source temperature using plume model
///
/// For the coarse-grid FFD,if the heat source cell is too large,the predicted
/// velocity at the heat source cell might not be accurate, the plume model can
/// be applied to predict that veolcity, and to correct the FFD prediction.
///
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param BINDEX Pointer to boundary index
///
///\return void No return needed
///////////////////////////////////////////////////////////////////////////////
void plume_thermal(PARA_DATA *para, REAL **var, int **BINDEX){
  int i, j, k;
  int it;
  int imax = para->geom->imax, jmax = para->geom->jmax;
  int kmax = para->geom->kmax;
  int indexmax,indexmax_us;
  int zone_num;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  REAL *x = var[X],  *y = var[Y],  *z = var[Z]; 
  REAL *gx = var[GX],  *gy = var[GY],  *gz = var[GZ]; 
  REAL *u = var[VX], *v = var[VY], *w = var[VZ];
  REAL *T=var[TEMP],*T0=var[TMP1];
  REAL *T_temp=var[TMP2];
  REAL *flagp = var[FLAGP];
  REAL Dx,Dy,Dz;
  REAL Prt=para->prob->Prt;
  REAL rho=para->prob->rho;
  REAL cp=para->prob->spec;
  REAL beta=para->prob->beta;
  REAL b,dT,dT1;
  REAL zv=para->geom->zv;
  REAL z_s;
  REAL Kv=0.122;
  REAL Kt=0.45; //0.4
  REAL Kq=0.005;
  REAL eta=0.09;
  REAL Pc;
  REAL T_amb;
  REAL mass=0;
  REAL area=0;
  REAL T_plume;

  indexmax=para->geom->index[1];
  indexmax_us=para->geom->index[2];


  for(it=indexmax+1;it<=indexmax_us;it++) {
    i=BINDEX[0][it];
    j=BINDEX[1][it];
    k=BINDEX[2][it];
    zone_num=BINDEX[16][FIX(i,j,k)]; 
    Pc=para->bc->t_bc[zone_num];
    BINDEX[20][FIX(i,j,k)]=1;
    Dx = gx[FIX(i,j,k)]-gx[FIX(i-1,j,k)];
    Dy = gy[FIX(i,j ,k)]-gy[FIX(i,j-1,k)];
    Dz = gz[FIX(i,j ,k)]-gz[FIX(i,j,k-1)];



    switch(para->prob->gravdir) {
      case GRAVX:
      case GRAVXN:
        z_s=x[FIX(i,j,k)]+zv;
        b=6/5*eta*(x[FIX(i,j,k)]-z_s);
        dT=Kt*pow(Pc,0.66667)*pow((x[FIX(i,j,k)]-z_s),-1.66667);
        dT1=b*b/(Dy*Dz)*(1-exp(-Dy*Dz/(b*b)))*dT;
        T_temp[FIX(i,j,k)]= T[FIX(i,j,k)];
        T_amb=0.25f*(T[FIX(i,j-1,k)]+T[FIX(i,j+1,k)]+T[FIX(i,j,k-1)]+T[FIX(i,j,k+1)]);
        T_plume=dT+T_amb;
        T0[FIX(i,j,k)] += (T_plume-T[FIX(i,j,k)]);

        break;

      case GRAVY:
      case GRAVYN:

        z_s=y[FIX(i,j,k)]+zv;
        b=6/5*eta*(y[FIX(i,j,k)]-z_s);
        dT=Kt*pow(Pc,0.66667)*pow((y[FIX(i,j,k)]-z_s),-1.66667);
        dT1=b*b/(Dx*Dz)*(1-exp(-Dx*Dz/(b*b)))*dT;
        T_temp[FIX(i,j,k)]= T[FIX(i,j,k)];
        T_amb=0.25f*(T[FIX(i-1,j,k)]+T[FIX(i+1,j,k)]+T[FIX(i,j,k-1)]+T[FIX(i,j,k+1)]);
        T_plume=dT+T_amb;
        T0[FIX(i,j,k)] += (T_plume-T[FIX(i,j,k)]);

        break;

      case GRAVZ:
      case GRAVZN:

        z_s=z[FIX(i,j,k)]+zv;
        b=6/5*eta*(z[FIX(i,j,k)]-z_s);
        dT=Kt*pow(Pc,0.66667)*pow((z[FIX(i,j,k)]-z_s),-1.66667);
        dT1=b*b/(Dx*Dy)*(1-exp(-Dx*Dy/(b*b)))*dT;
        T_temp[FIX(i,j,k)]= T[FIX(i,j,k)];
        T_amb=0.25f*(T[FIX(i-1,j,k)]+T[FIX(i+1,j,k)]+T[FIX(i,j-1,k)]+T[FIX(i,j+1,k)]);
        T_plume=dT+T_amb;
        T0[FIX(i,j,k)] += (T_plume-T[FIX(i,j,k)]);


        break;
    
    }
  }
}



