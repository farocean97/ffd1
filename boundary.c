#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "data_structure.h"
#include "boundary.h"
#include "utility.h"



void set_vel_bnd(PARA_DATA *para, REAL **var,int **BINDEX) {
  REAL *u = var[VX], *v = var[VY], *w = var[VZ];

  set_bnd_vel(para, var, VX, u, BINDEX);
  set_bnd_vel(para, var, VY, v, BINDEX);
  set_bnd_vel(para, var, VZ, w, BINDEX); 
}//End of set_vel_bnd() 


/******************************************************************************
|  Set the boundary conditions
******************************************************************************/
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


/**************************************************************************************************************
|  Set the boundary conditions for velocity
*******************************************************************************************************************/
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
    case VX:
      for(it=1;it<=indexmax;it++){
        i=BINDEX[0][it];
        j=BINDEX[1][it];
        k=BINDEX[2][it]; 
        if(flagp[FIX(i,j,k)]==0) continue;
        if(flagu[FIX(i,j,k)]==0){
          psi[FIX(i,j,k)]= para->bc->u_bc[BINDEX[13][FIX(i,j,k)]]; 
          b[FIX(i+1,j,k)] += para->bc->um_bc[BINDEX[13][FIX(i,j,k)]]*
                  (gy[FIX(i,j,k)]-gy[FIX(i,j-1,k)])*(gz[FIX(i,j,k)]-gz[FIX(i,j,k-1)]);
        }
        if(flagu[FIX(i,j,k)]>1) psi[FIX(i,j,k)]= 0;
        if(flagu[FIX(i,j,k)]==1 && i<=imax) {
          psi[FIX(i,j,k)]= psi[FIX(i+1,j,k)]; 
          ap[FIX(i+1,j,k)] -= aw[FIX(i+1,j,k)];
          aw[FIX(i+1,j,k)]=0; 
        }

        ii=max(i-1,0);
        if(flagu[FIX(ii,j,k)]==0){
          psi[FIX(ii,j,k)]= para->bc->u_bc[BINDEX[13][FIX(ii,j,k)]]; 
          b[FIX(ii-1,j,k)] += para->bc->um_bc[BINDEX[13][FIX(ii,j,k)]]*
                  (gy[FIX(i,j,k)]-gy[FIX(i,j-1,k)])*(gz[FIX(i,j,k)]-gz[FIX(i,j,k-1)]);
        }
        if(flagu[FIX(ii,j,k)]>1 ) psi[FIX(ii,j,k)]= 0;
        if(flagu[FIX(ii,j,k)]==1 && ii>0){
          psi[FIX(ii,j,k)]= psi[FIX(ii-1,j,k)]; 
          ap[FIX(ii-1,j,k)] -= ae[FIX(ii-1,j,k)]; 
          ae[FIX(ii-1,j,k)]=0; 
        }
      }

      for(it=1;it<=pindexmax_x;it++){
        i=BINDEX[17][it];
        j=BINDEX[18][it];
        k=BINDEX[19][it];
        psi[FIX(i,j,k)]=0;
      }
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

    case VY:
      for(it=1;it<=indexmax;it++) {
        i=BINDEX[0][it];
        j=BINDEX[1][it];
        k=BINDEX[2][it];
        if(flagp[FIX(i,j,k)]==0) continue; //heat source cell 
        if(flagv[FIX(i,j,k)]==0) psi[FIX(i,j,k)]= para->bc->v_bc[BINDEX[14][FIX(i,j,k)]]; 
        if(flagv[FIX(i,j,k)]> 1) psi[FIX(i,j,k)]= 0;
        if(flagv[FIX(i,j,k)]==1 && j<=jmax) {
          psi[FIX(i,j,k)]= psi[FIX(i,j+1,k)]; 
          ap[FIX(i,j+1,k)] -= as[FIX(i,j+1,k)]; 
          as[FIX(i,j+1,k)]=0; 
        }

        jj=max(j-1,0);
        if(flagv[FIX(i,jj,k)]==0) psi[FIX(i,jj,k)]= para->bc->v_bc[BINDEX[14][FIX(i,jj,k)]];
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
        
    case VZ:
      for(it=1;it<=indexmax;it++){
        i=BINDEX[0][it];
        j=BINDEX[1][it];
        k=BINDEX[2][it];
        if(flagp[FIX(i,j,k)]==0) continue;
        if(flagw[FIX(i,j,k)]==0) psi[FIX(i,j,k)]= para->bc->w_bc[BINDEX[15][FIX(i,j,k)]]; 
        if(flagw[FIX(i,j,k)]>1) psi[FIX(i,j,k)]= 0; 
        if(flagw[FIX(i,j,k)]==1 && k<=kmax){
          psi[FIX(i,j,k)]= psi[FIX(i,j,k+1)];
          ap[FIX(i,j,k+1)] -= ab[FIX(i,j,k+1)];
          ab[FIX(i,j,k+1)]=0; 
        }

        kk=max(k-1,0);
        if(flagw[FIX(i,j,kk)]==0) psi[FIX(i,j,kk)]= para->bc->w_bc[BINDEX[15][FIX(i,j,kk)]]; 
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
  int priority=0;
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
    
    /*******X direction EAST FACE********/
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
    
    /*******X direction WEST FACE********/
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
    
    /*******Y direction NORTH FACE********/
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
    
    
    /*******Y direction SOUTH FACE********/
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


    /*******Z direction FRONT FACE********/
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
    
    
    /*******Z direction BACK FACE********/
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
  }//end_for
  
  
  /**********************************Partitions****************************************************/
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
  
  /**********************************Heat sources****************************************************/
  for(it=indexmax+1;it<=indexmax_us;it++) {
    i=BINDEX[0][it];
    j=BINDEX[1][it];
    k=BINDEX[2][it];
    zone_num=BINDEX[16][FIX(i,j,k)]; 
    if(para->prob->plume_mod==0) {
      b[FIX(i,j,k)] += para->bc->t_bc[zone_num]/(rho*cp);
    }
    else {
      b[FIX(i,j,k)] += para->bc->t_bc[zone_num]/(rho*cp);

      //qflux= para->bc->qdiff;
      //af[FIX(i,j,k-1)]=0;
      //b[FIX(i,j,k-1)] += 1/(rho*cp)*qflux*(gy[FIX(i,j,k)]-gy[FIX(i,j-1,k)])
      //                                   *(gx[FIX(i,j,k)]-gx[FIX(i-1,j,k)]);
      //ab[FIX(i,j,k+1)]=0;
      //b[FIX(i,j,k+1)] += 1/(rho*cp)*qflux*(gy[FIX(i,j,k)]-gy[FIX(i,j-1,k)])
      //                                   *(gx[FIX(i,j,k)]-gx[FIX(i-1,j,k)]);
      //aw[FIX(i+1,j,k)]=0;
      //b[FIX(i+1,j,k)] += 1/(rho*cp)*qflux*(gy[FIX(i,j,k)]-gy[FIX(i,j-1,k)])
      //                                   *(gz[FIX(i,j,k+1)]-gz[FIX(i,j,k)]);
      //ae[FIX(i-1,j,k)]=0;
      //b[FIX(i-1,j,k)] += 1/(rho*cp)*qflux*(gy[FIX(i,j,k)]-gy[FIX(i,j-1,k)])
      //                                   *(gz[FIX(i,j,k+1)]-gz[FIX(i,j,k)]);
      //as[FIX(i,j+1,k)]=0;
      //b[FIX(i,j+1,k)] += 1/(rho*cp)*qflow_a;
  
    }
  }
  

}



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


/******************************************************************************
|  Set the boundary conditions for pressure
******************************************************************************/
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



/******************************************************************************
| Check the mass conservation of the domain
******************************************************************************/
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

  /*---------------------------------------------------------------------------
  | Compute the total inflow
  ---------------------------------------------------------------------------*/

  zone_inlet=para->geom->zone_inlet;
  zone_outlet=para->geom->zone_outlet;
  
  
  /***************inlet***************************************/
  
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
  
  
  
  
  /***************outlet***************************************/
  
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
  
  //printf("mass_in = %f, mass_out = %f, mass_in/mass_out=%f\n", mass_in, mass_out,mass_in/mass_out); 
  
  /***************mass correction***************************************/

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
          u[FIX(SI,j+1,k+1)] -= flagp[FIX(SI,j+1,k+1)]*(mass_in+mass_out)/area_out; 
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



void plume(PARA_DATA *para, REAL **var, int **BINDEX){
  int i, j, k;
  int ii,jj,kk;
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
  REAL alpha;
  REAL qflux_r,qflow_a;
  REAL b,dx,dz,dist,dT,radius;
  REAL a1,Dx,Dy,Dz;
  REAL vc,v_plume;
  REAL zv=para->geom->zv;
  REAL z_s;
  REAL Kv=0.122;
  REAL Kt=0.45; //0.4
  REAL Qp;
  REAL Kq=0.0045;
  REAL eta=0.09;
  REAL Pc;
  REAL T_amb;
  REAL T_s=300;
  REAL mass=0;
  REAL area=0;
  REAL *s_plume=var[ER];
  
  indexmax=para->geom->index[1];
  indexmax_us=para->geom->index[2];


  /**********************************Heat sources********************************/
  for(it=indexmax+1;it<=indexmax_us;it++) {
    i=BINDEX[0][it];
    j=BINDEX[1][it];
    k=BINDEX[2][it];
    zone_num=BINDEX[16][FIX(i,j,k)]; 

    para->geom->jplume=j;
    para->geom->iplume=i;
    para->geom->kplume=k;
    z_s=y[FIX(i,j,k)]+zv;

    Pc=para->bc->t_bc[zone_num];
    //printf("Pc=%f\n",Pc);
    
    //T_amb= 0.25f*(T[FIX(i+1,j,k)]+T[FIX(i-1,j,k)]+T[FIX(i,j,k+1)]+T[FIX(i,j,k-1)]);
    //T_amb= T[FIX(imax,j,k)];

    //zv= y[FIX(i,j,k)]-pow((T_s-T_amb),-0.6)*pow(Kt,0.6)*pow(Pc,0.4);
    //printf("zv=%f\n",zv);

    jj=0;

    for(jj=0;jj<para->geom->plmax;jj++) {
      if(j+jj>=jmax || flagv[FIX(i,j+jj,k)]>=0) continue;
      Dx = gx[FIX(i,j,k)]-gx[FIX(i-1,j,k)];
      Dy = gy[FIX(i,j ,k)]-gy[FIX(i,j-1,k)];
      Dz = gz[FIX(i,j ,k)]-gz[FIX(i,j,k-1)];
      b=6/5*eta*(gy[FIX(i,j+jj,k)]-z_s);
      a1=3.14f*b*b*(1-exp(-Dx*Dz/(b*b)))/(Dx*Dz);
      vc=Kv*pow(Pc,0.33333)*pow((gy[FIX(i,j+jj,k)]-z_s),-0.33333);
      Qp=Kq*pow(Pc,0.33333)*pow((gy[FIX(i,j+jj,k)]-z_s),1.666667);
      //v[FIX(i,j+jj,k)]=vc;
      //v[FIX(i,j+jj,k)]=a1*vc;//average velocity by massflow
      //v_plume=a1*vc;
      v_plume=0.63f*a1*vc;//0.63
      //v_plume=0.63f*Qp/(Dx*Dz);
      s_plume[FIX(i,j+jj,k)]=v_plume-v[FIX(i,j+jj,k)];

      //printf("vplume=%f\tV=%f\ts_plume=%f\n",Qp,v[FIX(i,j+jj,k)],s_plume[FIX(i,j+jj,k)]);

      //alpha=(beta+vt[FIX(i,j,k)])/Prt;
      //radius=gx[FIX(i,j,k)]-gx[FIX(i-1,j,k)];
      //qflux_r = 2.0f*radius/(b*b)*(T[FIX(i,j,k)]-T_amb)*(rho*cp*alpha)*exp(-radius*radius/(b*b));
      //qflow_a = 5.0f/3.0f*3.1415f*dT/(gy[FIX(i,j+jj,k)]-zv)*(1-exp(-radius*radius/(b*b)));
      //para->bc->qdiff =qflux_r;
      //para->bc->qflow_a=qflow_a;
      //para->bc->qflow_diff=qflow_a+2.0f*qflux_r*(gy[FIX(i,j,k)]-gy[FIX(i,j-1,k)])
      //                             *((gx[FIX(i,j,k)]-gx[FIX(i-1,j,k)])
      //                              +(gz[FIX(i,j,k)]-gz[FIX(i,j,k-1)]));
                                

    }


      
      /*if((j+jj)==jndex) {
        mass=0;
        area=0;
        for(ii=1;ii<=imax;ii++)
          for(kk=1;kk<=kmax;kk++) {
            // if(flagp[FIX(ii,j+jj,kk)]>=0) continue;
            dx=x[FIX(ii,j+jj,kk)]-x[FIX(i,j+jj,k)];
            dz=z[FIX(ii,j+jj,kk)]-z[FIX(i,j+jj,k)];
            dist=sqrt(dx*dx+dz*dz);
            if(dist<=b) {
              BINDEX[20][FIX(ii,j+jj,kk)]=1;
              v[FIX(ii,j+jj,kk)]=vc* exp(-dist*dist/(b*b));
              T[FIX(ii,j+jj,kk)]=dT* exp(-dist*dist/(b*b))+T_amb;
              mass += 3600*v[FIX(ii,j+jj,kk)]*(gx[FIX(ii,j+jj,kk)]-gx[FIX(ii-1,j+jj,kk)])
                                             *(gz[FIX(ii,j+jj,kk)]-gz[FIX(ii,j+jj,kk-1)]);
              area += (gx[FIX(ii,j+jj,kk)]-gx[FIX(ii-1,j+jj,kk)])
                     *(gz[FIX(ii,j+jj,kk)]-gz[FIX(ii,j+jj,kk-1)]);
            }
          }
      } 
    }*/
  }
  
  //printf("y=%f\tvc=%f\tmass=%f\tarea=%f\n",gy[FIX(1,jndex,1)],v[FIX(i_center,jndex
   //                                                          ,k_center)],mass,area );

}

void plume_thermal(PARA_DATA *para, REAL **var, int **BINDEX){
  int i, j, k;
  int ii,jj,kk;
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
  REAL alpha;
  REAL qflux_r,qflow_a;
  REAL b,dist,dT,dT1,radius;
  REAL vc;
 // REAL zv=para->geom->zv;
  REAL zv=-3.0f; //-6.0
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


  /**********************************Heat sources********************************/
  for(it=indexmax+1;it<=indexmax_us;it++) {
    i=BINDEX[0][it];
    j=BINDEX[1][it];
    k=BINDEX[2][it];
    zone_num=BINDEX[16][FIX(i,j,k)]; 
    z_s=y[FIX(i,j,k)]+zv;

    para->geom->jplume=j;
    para->geom->iplume=i;
    para->geom->kplume=k;
    Pc=para->bc->t_bc[zone_num];

    jj=0;
    
    //for(jj=0;jj<para->geom->plmax;jj++) {
    for(jj=0;jj<1;jj++) {
      BINDEX[20][FIX(i,j+jj,k)]=1;
      Dx = gx[FIX(i,j,k)]-gx[FIX(i-1,j,k)];
      Dy = gy[FIX(i,j ,k)]-gy[FIX(i,j-1,k)];
      Dz = gz[FIX(i,j ,k)]-gz[FIX(i,j,k-1)];
      b=6/5*eta*(y[FIX(i,j+jj,k)]-z_s);
      dT=Kt*pow(Pc,0.66667)*pow((y[FIX(i,j+jj,k)]-z_s),-1.66667);
      dT1=b*b/(Dx*Dz)*(1-exp(-Dx*Dz/(b*b)))*dT;
      T_temp[FIX(i,j+jj,k)]= T[FIX(i,j+jj,k)];
      //T_amb=T[FIX(imax,j+jj,k)];
      T_amb=T[FIX(i-1,j+jj,k)];
      //T[FIX(i,j+jj,k)]=dT+T_amb;
      T_plume=dT+T_amb;
      //T_plume=dT1+T_amb;
      T0[FIX(i,j+jj,k)] += (T_plume-T[FIX(i,j+jj,k)]);
      //T[FIX(i,j+jj,k)]=dT1+T_amb;
      //T[FIX(i,j+jj,k)]=299;
      //printf("Tplume=%f\t%f\t%f\n",T[FIX(i,j+jj,k)],T_plume, dT1);

     }
 

  }


}


