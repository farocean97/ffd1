///////////////////////////////////////////////////////////////////////////////
///
/// \file   utility.c
///
/// \brief  Various functions used for data processing in FFD
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
/// This file provides functions for processing the data in FFD
///
///////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "data_structure.h"
#include "utility.h"
#include "interpolation.h"


///////////////////////////////////////////////////////////////////////////////
///\brief  Conservation restoration of semi-Lagrangian
///
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param psi0 Pointer to the computed variables before advection
///\param psi Pointer to the computed variables after advection
///\param BINDEX Pointer to boundary index
///
///\return void No return needed
///////////////////////////////////////////////////////////////////////////////
void psi_conservation(PARA_DATA *para, REAL **var, REAL *psi,
                      REAL *psi0,int **BINDEX) {
  int i,j,k;
  int imax = para->geom->imax, jmax = para->geom->jmax, kmax=para->geom->kmax;
  REAL *u = var[VX], *v = var[VY];
  REAL *gx = var[GX], *gy = var[GY], *gz = var[GZ];
  REAL *flagp = var[FLAGP];
  REAL dA;
  REAL dt=para->mytime->dt;
  REAL mass0=0, mass=0.0000001f,massstar=0;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  REAL qin,qout;
  REAL area=0;
  REAL *dens_s=var[DENS];
  REAL dens=0.0000001f, dens1=0.0000001f;
  REAL eta;
  REAL rho=para->prob->rho;
  REAL cp=para->prob->spec;

  if(para->prob->plume_mod==1) {
    mass0=0;
    FOR_EACH_CELL
      if(flagp[FIX(i,j,k)]<=0) {
        var[LOCMIN][FIX(i,j,k)]=check_min(para, psi0, BINDEX[4][FIX(i,j,k)], 
                                BINDEX[5][FIX(i,j,k)],  BINDEX[6][FIX(i,j,k)]); 
        var[LOCMAX][FIX(i,j,k)]=check_max(para, psi0, BINDEX[4][FIX(i,j,k)], 
                                BINDEX[5][FIX(i,j,k)],  BINDEX[6][FIX(i,j,k)]); 
      }
    END_FOR

    qin= inflow(para,var,psi0,BINDEX);
    qout=outflow(para,var,psi0,BINDEX);
    mass0 +=(qin+qout)*dt;

    FOR_EACH_CELL
      if(flagp[FIX(i,j,k)]>0) continue;
      dA= (gx[FIX(i,j,k)]-gx[FIX(i-1,j,k)])*(gz[FIX(i,j,k)]-gz[FIX(i,j,k-1)])
          *(gy[FIX(i,j,k)]-gy[FIX(i,j-1,k)]);
      area += dA;
      if(BINDEX[20][FIX(i,j,k)]==1) mass0 += rho*cp*var[TMP2][FIX(i,j,k)]*dA;
      else mass0 += rho*cp*psi0[FIX(i,j,k)]*dA;
      massstar += rho*cp*psi0[FIX(i,j,k)]*dA;
      mass  += rho*cp*psi[FIX(i,j,k)]*dA;
      dens +=  (float) fabs(psi[FIX(i,j,k)]-var[LOCMIN][FIX(i,j,k)]);
      dens1 +=  (float) fabs(psi[FIX(i,j,k)]-var[LOCMAX][FIX(i,j,k)]);   
    END_FOR

    eta=mass0-mass;
    if(eta<0) {   // mass generation occured
      FOR_EACH_CELL
        if(flagp[FIX(i,j,k)]>0 ) continue;
        dens_s[FIX(i,j,k)]=(psi[FIX(i,j,k)]-var[LOCMIN][FIX(i,j,k)])/dens
                       *eta/(rho*cp*(gx[FIX(i,j,k)]-gx[FIX(i-1,j,k)])
                                   *(gz[FIX(i,j,k)]-gz[FIX(i,j,k-1)])
                                   *(gy[FIX(i,j,k)]-gy[FIX(i,j-1,k)]));
        psi[FIX(i,j,k)] += (float) fabs(psi[FIX(i,j,k)]-var[LOCMIN][FIX(i,j,k)])/dens
                       *eta/(rho*cp*(gx[FIX(i,j,k)]-gx[FIX(i-1,j,k)])
                                   *(gz[FIX(i,j,k)]-gz[FIX(i,j,k-1)])
                                   *(gy[FIX(i,j,k)]-gy[FIX(i,j-1,k)]));
      END_FOR
    }
    else {
      FOR_EACH_CELL
        if(flagp[FIX(i,j,k)]>0 ) continue;
          dens_s[FIX(i,j,k)]=(psi[FIX(i,j,k)]-var[LOCMAX][FIX(i,j,k)])/dens1
                       *eta/(rho*cp*(gx[FIX(i,j,k)]-gx[FIX(i-1,j,k)])
                                   *(gz[FIX(i,j,k)]-gz[FIX(i,j,k-1)])
                                   *(gy[FIX(i,j,k)]-gy[FIX(i,j-1,k)]));
          psi[FIX(i,j,k)] += (float) fabs(psi[FIX(i,j,k)]-var[LOCMAX][FIX(i,j,k)])/dens1
                      *eta/(rho*cp*(gx[FIX(i,j,k)]-gx[FIX(i-1,j,k)])
                                  *(gz[FIX(i,j,k)]-gz[FIX(i,j,k-1)])
                                  *(gy[FIX(i,j,k)]-gy[FIX(i,j-1,k)]));
      END_FOR
    }
  }
  else { // without plume model
    mass0=0;
    FOR_EACH_CELL
      if(flagp[FIX(i,j,k)]<=0)  {
      var[LOCMIN][FIX(i,j,k)]=check_min(para, psi0, BINDEX[4][FIX(i,j,k)], 
                                      BINDEX[5][FIX(i,j,k)],  BINDEX[6][FIX(i,j,k)]); 
      var[LOCMAX][FIX(i,j,k)]=check_max(para, psi0, BINDEX[4][FIX(i,j,k)], 
                                      BINDEX[5][FIX(i,j,k)],  BINDEX[6][FIX(i,j,k)]); 
      }
    END_FOR

    qin= inflow(para,var,psi0,BINDEX);
    qout=outflow(para,var,psi0,BINDEX);
    mass0 +=(qin+qout)*dt;

    FOR_EACH_CELL
      if(flagp[FIX(i,j,k)]>0) continue;
      dA= (gx[FIX(i,j,k)]-gx[FIX(i-1,j,k)])*(gz[FIX(i,j,k)]-gz[FIX(i,j,k-1)])
                                           *(gy[FIX(i,j,k)]-gy[FIX(i,j-1,k)]);
      area += dA;
      mass0 += rho*cp*psi0[FIX(i,j,k)]*dA;
      massstar += rho*cp*psi0[FIX(i,j,k)]*dA;
      mass  += rho*cp*psi[FIX(i,j,k)]*dA;
      dens +=  (float) fabs(psi[FIX(i,j,k)]-var[LOCMIN][FIX(i,j,k)]);
      dens1 +=  (float) fabs(psi[FIX(i,j,k)]-var[LOCMAX][FIX(i,j,k)]);   
    END_FOR

    massstar =mass0;
    eta=mass0-mass;
    if(eta<0) {   // mass generation occured
      FOR_EACH_CELL
        if(flagp[FIX(i,j,k)]>0) continue;
        dens_s[FIX(i,j,k)]=(psi[FIX(i,j,k)]-var[LOCMIN][FIX(i,j,k)])/dens
                        *eta/(rho*cp*(gx[FIX(i,j,k)]-gx[FIX(i-1,j,k)])
                                    *(gz[FIX(i,j,k)]-gz[FIX(i,j,k-1)])
                                    *(gy[FIX(i,j,k)]-gy[FIX(i,j-1,k)]));
        psi[FIX(i,j,k)] += (float) fabs(psi[FIX(i,j,k)]-var[LOCMIN][FIX(i,j,k)])/dens
                        *eta/(rho*cp*(gx[FIX(i,j,k)]-gx[FIX(i-1,j,k)])
                                    *(gz[FIX(i,j,k)]-gz[FIX(i,j,k-1)])
                                    *(gy[FIX(i,j,k)]-gy[FIX(i,j-1,k)]));
      END_FOR
    }
    else {
      FOR_EACH_CELL
        if(flagp[FIX(i,j,k)]>0) continue;
          dens_s[FIX(i,j,k)]=(psi[FIX(i,j,k)]-var[LOCMAX][FIX(i,j,k)])/dens1
                        *eta/(rho*cp*(gx[FIX(i,j,k)]-gx[FIX(i-1,j,k)])
                                    *(gz[FIX(i,j,k)]-gz[FIX(i,j,k-1)])
                                    *(gy[FIX(i,j,k)]-gy[FIX(i,j-1,k)]));
          psi[FIX(i,j,k)] += (float) fabs(psi[FIX(i,j,k)]-var[LOCMAX][FIX(i,j,k)])/dens1
                        *eta/(rho*cp*(gx[FIX(i,j,k)]-gx[FIX(i-1,j,k)])
                                    *(gz[FIX(i,j,k)]-gz[FIX(i,j,k-1)])
                                    *(gy[FIX(i,j,k)]-gy[FIX(i,j-1,k)]));
      END_FOR
    }
  } 


} // End of psi_conservation()



///////////////////////////////////////////////////////////////////////////////
///\brief  Compute the energy or specices flow out of the domain
///
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param psi Pointer to the computed variables before advection
///\param BINDEX Pointer to boundary index
///
///\return outflow rate
///////////////////////////////////////////////////////////////////////////////
REAL outflow(PARA_DATA *para, REAL **var, REAL *psi, int **BINDEX)
{
  int i, j, k;
  int imax = para->geom->imax, jmax = para->geom->jmax;
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  REAL *gx = var[GX], *gy = var[GY], *gz = var[GZ];
  REAL *u = var[VX], *v = var[VY], *w = var[VZ];
  REAL mass_out=0;
  REAL *flagp = var[FLAGP];
  int zone_inlet,zone_outlet,zone_num;
  int SI,EI,SJ,EJ,SK,EK;
  REAL rho=para->prob->rho;
  REAL cp=para->prob->spec;

  zone_inlet=para->geom->zone_inlet;
  zone_outlet=para->geom->zone_outlet;

  for(zone_num=zone_inlet+1;zone_num<=zone_outlet;zone_num++) {
    SI=BINDEX[7][zone_num];
    EI=BINDEX[8][zone_num];
    SJ=BINDEX[9][zone_num];
    EJ=BINDEX[10][zone_num];
    SK=BINDEX[11][zone_num];
    EK=BINDEX[12][zone_num];

    if(SI==EI){
      for(j=SJ ;j<EJ ;j++)
        for(k=SK ;k<EK ;k++){
          mass_out += rho*cp*psi[FIX(SI,j+1,k+1)]*flagp[FIX(SI,j+1,k+1)]*u[FIX(SI,j+1,k+1)]
                     *(gy[FIX(SI,j+1,k+1)]-gy[FIX(SI,j,k+1)])* (gz[FIX(SI,j+1,k+1)]-gz[FIX(SI,j+1,k)]);
        }

    }
    if(SJ==EJ){
      for(i=SI ;i<EI ;i++)
        for(k=SK ;k<EK ;k++) {
          mass_out += rho*cp*psi[FIX(i+1,SJ,k+1)]*flagp[FIX(i+1,SJ,k+1)]*v[FIX(i+1,SJ,k+1)]
                     *(gx[FIX(i+1,SJ,k+1)]-gx[FIX(i,SJ,k+1)])* (gz[FIX(i+1,SJ,k+1)]-gz[FIX(i+1,SJ,k)]);
        }

    }

    if(SK==EK){
      for(i=SI ;i<EI ;i++)
        for(j=SJ ;j<EJ ;j++){
          mass_out += rho*cp*psi[FIX(i+1,j+1,SK)]*flagp[FIX(i+1,j+1,SK)]*w[FIX(i+1,j+1,SK)]
                    *(gx[FIX(i+1,j+1,SK)]-gx[FIX(i,j+1,SK)])* (gy[FIX(i+1,j+1,SK)]-gy[FIX(i+1,j,SK)]);
        }
    }
  }

  return mass_out;

}// End of outflow()



///////////////////////////////////////////////////////////////////////////////
///\brief  Compute the energy or specices flow into the domain
///
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param psi Pointer to the computed variables before advection
///\param BINDEX Pointer to boundary index
///
///\return inflow rate
///////////////////////////////////////////////////////////////////////////////
REAL inflow(PARA_DATA *para, REAL **var, REAL *psi, int **BINDEX){
  int i, j, k;
  int imax = para->geom->imax, jmax = para->geom->jmax;
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  REAL *gx = var[GX], *gy = var[GY], *gz = var[GZ];
  REAL *u = var[VX], *v = var[VY], *w = var[VZ];
  REAL mass_in=0;
  REAL *flagp = var[FLAGP];
  int zone_inlet,zone_outlet,zone_num;
  int SI,EI,SJ,EJ,SK,EK;
  REAL rho=para->prob->rho;
  REAL cp=para->prob->spec;


  zone_inlet=para->geom->zone_inlet;
  zone_outlet=para->geom->zone_outlet;


  for(zone_num=1;zone_num<=zone_inlet;zone_num++){
    SI=BINDEX[7][zone_num];
    EI=BINDEX[8][zone_num];
    SJ=BINDEX[9][zone_num];
    EJ=BINDEX[10][zone_num];
    SK=BINDEX[11][zone_num];
    EK=BINDEX[12][zone_num];

    if(SI==EI){
      for(j=SJ ;j<EJ ;j++)
        for(k=SK ;k<EK ;k++){
          mass_in +=rho*cp*para->bc->t_bc[zone_num]*flagp[FIX(SI,j+1,k+1)]*u[FIX(SI,j+1,k+1)]
                  *(gy[FIX(SI,j+1,k+1)]-gy[FIX(SI,j,k+1)])* (gz[FIX(SI,j+1,k+1)]-gz[FIX(SI,j+1,k)]);
        }
    }

    if(SJ==EJ){
      for(i=SI ;i<EI ;i++)
        for(k=SK ;k<EK ;k++){
          mass_in +=rho*cp*para->bc->t_bc[zone_num]* flagp[FIX(i+1,SJ,k+1)]*v[FIX(i+1,SJ,k+1)]
                  *(gx[FIX(i+1,SJ,k+1)]-gx[FIX(i,SJ,k+1)])* (gz[FIX(i+1,SJ,k+1)]-gz[FIX(i+1,SJ,k)]);
        }
    }

    if(SK==EK){
      for(i=SI ;i<EI ;i++)
        for(j=SJ ;j<EJ ;j++){
           mass_in += rho*cp*para->bc->t_bc[zone_num]*flagp[FIX(i+1,j+1,SK)]*w[FIX(i+1,j+1,SK)]
                     *(gx[FIX(i+1,j+1,SK)]-gx[FIX(i,j+1,SK)])* (gy[FIX(i+1,j+1,SK)]-gy[FIX(i+1,j,SK)]);
        }
    }
  }
  return mass_in;
} // End of inflow()

///////////////////////////////////////////////////////////////////////////////
///\brief  Compute the minimum value around the departure point
///
///
///\param para Pointer to FFD parameters
///\param psi Pointer to the computed variables before advection
///\param ci I-index of the departure cell
///\param cj J-index of the departure cell
///\param ck K-index of the departure cell
///
///\return minimum value
///////////////////////////////////////////////////////////////////////////////
REAL check_min(PARA_DATA *para, REAL *psi, int ci,int cj,int ck){
  int imax = para->geom->imax, jmax = para->geom->jmax; 
  int kmax = para->geom->kmax;
  int i, j, k;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  REAL tmp=psi[FIX(ci,cj,ck)];

  for(i=0;i<=1;i++)
    for(j=0;j<=1;j++)
          for(k=0;k<=1;k++)
        {
          if(tmp>psi[FIX(ci+i,cj+j,ck+k)]) tmp=psi[FIX(ci+i,cj+j,ck+k)];
        
        }

 return tmp;

}// End of check_min( )

///////////////////////////////////////////////////////////////////////////////
///\brief  Compute the maximum value around the departure point
///
///
///\param para Pointer to FFD parameters
///\param psi Pointer to the computed variables before advection
///\param ci I-index of the departure cell
///\param cj J-index of the departure cell
///\param ck K-index of the departure cell
///
///\return maximum value
///////////////////////////////////////////////////////////////////////////////
REAL check_max( PARA_DATA *para, REAL *psi, int ci,int cj,int ck) {
  int imax = para->geom->imax, jmax = para->geom->jmax; 
  int kmax = para->geom->kmax;
  int i, j, k;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  REAL tmp=psi[FIX(ci,cj,ck)];

  for(i=0;i<=1;i++)
    for(j=0;j<=1;j++)
          for(k=0;k<=1;k++)
        {
          if(tmp<psi[FIX(ci+i,cj+j,ck+k)]) tmp=psi[FIX(ci+i,cj+j,ck+k)];
        
        }
    
  return tmp;

}// End of check_max( )


///////////////////////////////////////////////////////////////////////////////
///\brief  Compute the minimum value of eight adjacent cells
///
///
///\param para Pointer to FFD parameters
///\param psi Pointer to the variable
///
///\return minimum value
///////////////////////////////////////////////////////////////////////////////
REAL check_min_pix(PARA_DATA *para, REAL *psi){
  int i, j, k;
  REAL tmp=psi[PIX(0,0,0)];

  for(i=0;i<=1;i++)
    for(j=0;j<=1;j++)
          for(k=0;k<=1;k++)
        {
          if(tmp>psi[PIX(i,j,k)]) tmp=psi[PIX(i,j,k)];
        
        }
  return tmp;

}// End of check_min_pix( )

///////////////////////////////////////////////////////////////////////////////
///\brief  Compute the maximum value of eight adjacent cells
///
///
///\param para Pointer to FFD parameters
///\param psi Pointer to the variable
///
///\return maximum value
///////////////////////////////////////////////////////////////////////////////
REAL check_max_pix( PARA_DATA *para, REAL *psi) {
  int i, j, k;
  REAL tmp=psi[PIX(0,0,0)];

  for(i=0;i<=1;i++)
    for(j=0;j<=1;j++)
       for(k=0;k<=1;k++) {
          if(tmp<psi[PIX(i,j,k)]) tmp=psi[PIX(i,j,k)];
        
        }
  return tmp;

}// End of check_max_pix()


///////////////////////////////////////////////////////////////////////////////
///\brief  Compute the minimum value of eight adjacent cells
///
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param BINDEX Pointer to boundary index
///
///\return energy through the wall
///////////////////////////////////////////////////////////////////////////////
REAL qwall(PARA_DATA *para, REAL **var,int **BINDEX) {
  int i, j, k;
  int it;
  int indexmax,indexmax_us,zone_num;
  int imax = para->geom->imax, jmax = para->geom->jmax;
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  REAL *psi=var[TEMP];
  REAL *gx = var[GX], *gy = var[GY], *gz = var[GZ];
  REAL coeff_h=para->prob->coeff_h;
  REAL rho=para->prob->rho;
  REAL cp=para->prob->spec;
  REAL qwall=0,qwall_temp=0;
  REAL qw=0,qe=0,qn=0,qs=0,qf=0,qb=0;
  REAL qflux;
  REAL *flagp = var[FLAGP],*flagu = var[FLAGU],*flagv = var[FLAGV],*flagw = var[FLAGW];

  indexmax=para->geom->index[1];
  indexmax_us=para->geom->index[2];

  for(it=1;it<=indexmax;it++){
    i=BINDEX[0][it];
    j=BINDEX[1][it];
    k=BINDEX[2][it];
    
    /*******X direction EAST FACE********/

    if(i<imax && flagp[FIX(i+1,j,k)]<0) {
      switch((int)flagu[FIX(i,j,k)]) {
        case 0: //inlet_surface
        case 1: //outlet_surface
          qwall_temp=0;
          break;
        case 2: //wall_surface
          zone_num=BINDEX[13][FIX(i,j,k)];
          if(para->bc->fltmp[zone_num]==0){
            qflux=para->bc->t_bc[zone_num];
            qwall_temp = qflux*(gy[FIX(i,j,k)]-gy[FIX(i,j-1,k)])*(gz[FIX(i,j,k)]-gz[FIX(i,j,k-1)]);
            qwall += qwall_temp;
          }
          else{
            qwall_temp=rho*cp*(para->bc->t_bc[zone_num]-psi[FIX(i+1,j,k)])*coeff_h
                              *(gy[FIX(i,j,k)]-gy[FIX(i,j-1,k)])*(gz[FIX(i,j,k)]-gz[FIX(i,j,k-1)]);
            qwall += qwall_temp; 
          }
          break;
  
        case 3://block_wall_surface
          zone_num=BINDEX[16][FIX(i,j,k)];
          if(para->bc->fltmp[zone_num]==0)
               {
               qflux=para->bc->t_bc[zone_num];
               qwall_temp = qflux*(gy[FIX(i,j,k)]-gy[FIX(i,j-1,k)])*(gz[FIX(i,j,k)]-gz[FIX(i,j,k-1)]);
               qwall += qwall_temp;
         
               }
              else
               {
               qwall_temp = rho*cp*(para->bc->t_bc[zone_num]-psi[FIX(i+1,j,k)])*coeff_h
                             *(gy[FIX(i,j,k)]-gy[FIX(i,j-1,k)])*(gz[FIX(i,j,k)]-gz[FIX(i,j,k-1)]);    
               qwall += qwall_temp;
               }
               break;
           }
           qe +=qwall_temp;
           }
            


      /*******X direction WEST FACE********/
        if(i>0 && flagp[FIX(i-1,j,k)]<0)
           {
           switch( (int)flagu[FIX(i-1,j,k)])
           {
              case 0:
            case 1:
              qwall_temp=0;
             break;
            case 2:
                         zone_num=BINDEX[13][FIX(i-1,j,k)];
              if(para->bc->fltmp[zone_num]==0)
               {
               qflux=para->bc->t_bc[zone_num];
               qwall_temp = qflux*(gy[FIX(i,j,k)]-gy[FIX(i,j-1,k)])*(gz[FIX(i,j,k)]-gz[FIX(i,j,k-1)]);
               qwall += qwall_temp;
       
               }
              else
               {
               qwall_temp = rho*cp*(para->bc->t_bc[zone_num]-psi[FIX(i-1,j,k)])*coeff_h
                           *(gy[FIX(i,j,k)]-gy[FIX(i,j-1,k)])*(gz[FIX(i,j,k)]-gz[FIX(i,j,k-1)]);
                      
               qwall += qwall_temp;                  
              }
              break;
            case 3:
                           zone_num=BINDEX[16][FIX(i,j,k)];
               if(para->bc->fltmp[zone_num]==0)
               {
               qflux=para->bc->t_bc[zone_num];
               qwall_temp = qflux*(gy[FIX(i,j,k)]-gy[FIX(i,j-1,k)])*(gz[FIX(i,j,k)]-gz[FIX(i,j,k-1)]);
               qwall += qwall_temp;         
               }
              else
               {
               qwall_temp = rho*cp*(para->bc->t_bc[zone_num]-psi[FIX(i-1,j,k)])*coeff_h
                          *(gy[FIX(i,j,k)]-gy[FIX(i,j-1,k)])*(gz[FIX(i,j,k)]-gz[FIX(i,j,k-1)]);
             
               qwall += qwall_temp;
               }
               break;
           }
           qw +=qwall_temp;
           }


      /*******Y direction NORTH FACE********/

        if(j<jmax && flagp[FIX(i,j+1,k)]<0)
           {
           switch( (int)flagv[FIX(i,j,k)])
           {
              case 0:
            case 1:
              qwall_temp=0;
             break;
            case 2:
                         zone_num=BINDEX[14][FIX(i,j,k)];
              if(para->bc->fltmp[zone_num]==0)
               {
               qflux=para->bc->t_bc[zone_num];
               qwall_temp = qflux*(gx[FIX(i,j,k)]-gx[FIX(i-1,j,k)])*(gz[FIX(i,j,k)]-gz[FIX(i,j,k-1)]);
               qwall += qwall_temp;       
               }
              else
               {
               qwall_temp = rho*cp*(para->bc->t_bc[zone_num]-psi[FIX(i,j+1,k)])*coeff_h
                         *(gx[FIX(i,j,k)]-gx[FIX(i-1,j,k)])*(gz[FIX(i,j,k)]-gz[FIX(i,j,k-1)]);
               qwall += qwall_temp;           
               }
              break;
            case 3:
                           zone_num=BINDEX[16][FIX(i,j,k)];
               if(para->bc->fltmp[zone_num]==0)
               {
               qflux=para->bc->t_bc[zone_num];
               qwall_temp = qflux*(gx[FIX(i,j,k)]-gx[FIX(i-1,j,k)])*(gz[FIX(i,j,k)]-gz[FIX(i,j,k-1)]);
               qwall += qwall_temp;
               }
              else
               {
               qwall_temp = rho*cp*(para->bc->t_bc[zone_num]-psi[FIX(i,j+1,k)])*coeff_h
                           *(gx[FIX(i,j,k)]-gx[FIX(i-1,j,k)])*(gz[FIX(i,j,k)]-gz[FIX(i,j,k-1)]);

               qwall += qwall_temp;
          
               }
               break;
           }
           qn +=qwall_temp;
           }


      /*******Y direction SOUTH FACE********/

        if(j>0 && flagp[FIX(i,j-1,k)]<0)
           {

           switch( (int)flagv[FIX(i,j-1,k)])
           {
            case 0:
            case 1:
              qwall_temp=0;
             break;
            case 2:
                         zone_num=BINDEX[14][FIX(i,j-1,k)];
              if(para->bc->fltmp[zone_num]==0)
               {
               qflux=para->bc->t_bc[zone_num];
               qwall_temp = qflux*(gx[FIX(i,j,k)]-gx[FIX(i-1,j,k)])*(gz[FIX(i,j,k)]-gz[FIX(i,j,k-1)]);    
               qwall += qwall_temp;
               }
              else
               {
                 qwall_temp = rho*cp*(para->bc->t_bc[zone_num]-psi[FIX(i,j-1,k)])*coeff_h
                        *(gx[FIX(i,j,k)]-gx[FIX(i-1,j,k)])*(gz[FIX(i,j,k)]-gz[FIX(i,j,k-1)]);
                 qwall += qwall_temp;           
               }
              break;
            case 3:
               zone_num=BINDEX[16][FIX(i,j,k)];
               if(para->bc->fltmp[zone_num]==0)
               {
               qflux=para->bc->t_bc[zone_num];
               qwall_temp = qflux*(gx[FIX(i,j,k)]-gx[FIX(i-1,j,k)])*(gz[FIX(i,j,k)]-gz[FIX(i,j,k-1)]); 
               qwall += qwall_temp;
               }
              else
               {
                qwall_temp = rho*cp*(para->bc->t_bc[zone_num]-psi[FIX(i,j-1,k)])*coeff_h
                             *(gx[FIX(i,j,k)]-gx[FIX(i-1,j,k)])*(gz[FIX(i,j,k)]-gz[FIX(i,j,k-1)]);
               qwall += qwall_temp;
          
               }
               break;
           }
           qs +=qwall_temp;
           }

      /*******Z direction FRONT FACE********/

        if(k<kmax && flagp[FIX(i,j,k+1)]<0)
           {
           switch( (int)flagw[FIX(i,j,k)])
           {
              case 0:
            case 1:
              qwall_temp=0;
             break;
            case 2:
                         zone_num=BINDEX[15][FIX(i,j,k)];
              if(para->bc->fltmp[zone_num]==0)
               {
               qflux=para->bc->t_bc[zone_num];
               qwall_temp = qflux*(gy[FIX(i,j,k)]-gy[FIX(i,j-1,k)])*(gx[FIX(i,j,k)]-gx[FIX(i-1,j,k)]);
               qwall += qwall_temp;
               }
              else
               {
                qwall_temp = rho*cp*(para->bc->t_bc[zone_num]-psi[FIX(i,j,k+1)])*coeff_h
                             *(gy[FIX(i,j,k)]-gy[FIX(i,j-1,k)])*(gx[FIX(i,j,k)]-gx[FIX(i-1,j,k)]);
                      
                qwall += qwall_temp;           
               }
              break;
            case 3:
                           zone_num=BINDEX[16][FIX(i,j,k)];
               if(para->bc->fltmp[zone_num]==0)
               {
               qflux=para->bc->t_bc[zone_num];
               qwall_temp = qflux*(gy[FIX(i,j,k)]-gy[FIX(i,j-1,k)])*(gx[FIX(i,j,k)]-gx[FIX(i-1,j,k)]);
               qwall += qwall_temp;
               }
              else
               {
                qwall_temp = rho*cp*(para->bc->t_bc[zone_num]-psi[FIX(i,j,k+1)])*coeff_h
                                *(gy[FIX(i,j,k)]-gy[FIX(i,j-1,k)])*(gx[FIX(i,j,k)]-gx[FIX(i-1,j,k)]);
           
               qwall += qwall_temp;
               }
               break;
           }
           qf +=qwall_temp;
           }

      /*******Z direction BACK FACE********/

        if(k>0 && flagp[FIX(i,j,k-1)]<0)
           {
           switch( (int)flagw[FIX(i,j,k-1)])
           {
              case 0:
            case 1:
              qwall_temp=0;
             break;
            case 2:
                         zone_num=BINDEX[15][FIX(i,j,k-1)];
              if(para->bc->fltmp[zone_num]==0)
               {
               qflux=para->bc->t_bc[zone_num];
               qwall_temp = qflux*(gy[FIX(i,j,k)]-gy[FIX(i,j-1,k)])*(gx[FIX(i,j,k)]-gx[FIX(i-1,j,k)]);  
               qwall += qwall_temp;
               }
              else
               {
               qwall_temp = rho*cp*(para->bc->t_bc[zone_num]-psi[FIX(i,j,k-1)])*coeff_h
                         *(gy[FIX(i,j,k)]-gy[FIX(i,j-1,k)])*(gx[FIX(i,j,k)]-gx[FIX(i-1,j,k)]);
            
               qwall += qwall_temp;           
               }
              break;
            case 3:
                           zone_num=BINDEX[16][FIX(i,j,k)];
               if(para->bc->fltmp[zone_num]==0)
               {
               qflux=para->bc->t_bc[zone_num];
               qwall_temp = qflux*(gy[FIX(i,j,k)]-gy[FIX(i,j-1,k)])*(gx[FIX(i,j,k)]-gx[FIX(i-1,j,k)]);     
               qwall += qwall_temp;
               }
              else
               {
                qwall_temp = rho*cp*(para->bc->t_bc[zone_num]-psi[FIX(i,j,k-1)])*coeff_h
                           *(gy[FIX(i,j,k)]-gy[FIX(i,j-1,k)])*(gx[FIX(i,j,k)]-gx[FIX(i-1,j,k)]);
            
               qwall += qwall_temp;
          
               }
               break;
           }
           qb +=qwall_temp;
      }

    }

  qwall += para->bc->qs;   //para->bc->t_bc[zone_num];

  return qwall;
      

} // End of qwall()


///////////////////////////////////////////////////////////////////////////////
///\brief  Swap the stored velocities
///
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///
///\return void No return needed
///////////////////////////////////////////////////////////////////////////////
void swapuvw(PARA_DATA *para, REAL **var) {
  int i, j, k;
  int imax = para->geom->imax, jmax = para->geom->jmax;
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
   REAL *u1 = var[VXM], *v1 = var[VYM], *w1 = var[VZM];  
  REAL *u = var[VX], *v = var[VY], *w = var[VZ];  
  

  FOR_ALL_CELL
    u[FIX(i,j,k)]=u1[FIX(i,j,k)];
    v[FIX(i,j,k)]=v1[FIX(i,j,k)];
    w[FIX(i,j,k)]=w1[FIX(i,j,k)];
  END_FOR

}// End of swapuvw

///////////////////////////////////////////////////////////////////////////////
///\brief  Compute the CFL number in the domain
///
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///
///\return void No return needed
///////////////////////////////////////////////////////////////////////////////
void cfl(PARA_DATA *para, REAL **var) {
  int i, j, k;
  int imax = para->geom->imax, jmax = para->geom->jmax;
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  REAL dt=para->mytime->dt;
  REAL *u = var[VX], *v = var[VY], *w = var[VZ];  
  REAL *cflx = var[CFLX], *cfly = var[CFLY],*cflz = var[CFLZ]; 
  REAL *gx = var[GX], *gy = var[GY], *gz = var[GZ];  

  FOR_EACH_CELL
  cflx[FIX(i,j,k)]=fabs(u[FIX(i,j,k)])*dt/(gx[FIX(i,j,k)]-gx[FIX(i-1,j,k)]);
  cfly[FIX(i,j,k)]=fabs(v[FIX(i,j,k)])*dt/(gy[FIX(i,j,k)]-gy[FIX(i,j-1,k)]);
  cflz[FIX(i,j,k)]=fabs(w[FIX(i,j,k)])*dt/(gz[FIX(i,j,k)]-gz[FIX(i,j,k-1)]);
  END_FOR

}// End of cfl()

///////////////////////////////////////////////////////////////////////////////
///\brief  Interpolate U from surface cell to grid point
///
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param u   Pointer to U at grid point
///\param u0   Pointer to U at surface cell
///
///\return void No return needed
///////////////////////////////////////////////////////////////////////////////
void u_vortex(PARA_DATA *para, REAL **var, REAL *u, REAL *u0)
{
  int i, j, k;
  int imax = para->geom->imax, jmax = para->geom->jmax;
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);  
  REAL *gx = var[GX], *gy = var[GY], *gz = var[GZ];  
  REAL *x = var[X], *y = var[Y], *z = var[Z];  
  REAL y_1,z_1;


    for(i=0;i<=imax;i++)
          for(j=0;j<=jmax;j++)
       for(k=0;k<=kmax;k++)
        {
          if(j==0 || j==jmax) u[FIX(i,j,k)]= (u0[FIX(i,j,k)]+u0[FIX(i,j,k+1)])/2.0f;
        else if (k==0 || k==kmax) u[FIX(i,j,k)]= (u0[FIX(i,j,k)]+u0[FIX(i,j+1,k)])/2.0f;
        else 
        {
           y_1= (gy[FIX(i,j,k)]-y[FIX(i,j,k)])/(y[FIX(i,j+1,k)]-y[FIX(i,j,k)]);
                   z_1= (gz[FIX(i,j,k)]-z[FIX(i,j,k)])/(z[FIX(i,j,k+1)]-z[FIX(i,j,k)]);
           u[FIX(i,j,k)]= interpolation_linear(y_1, z_1, u0[FIX(i,j,k)], u0[FIX(i,j,k+1)],
                                                         u0[FIX(i,j+1,k)], u0[FIX(i,j+1,k+1)]);
        
        }
        }

}// End of u_vortex()

///////////////////////////////////////////////////////////////////////////////
///\brief  Interpolate V from surface cell to grid point
///
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param v   Pointer to V at grid point
///\param v0   Pointer to V at surface cell
///
///\return void No return needed
///////////////////////////////////////////////////////////////////////////////
void v_vortex(PARA_DATA *para, REAL **var, REAL *v, REAL *v0)
{
  int i, j, k;
  int imax = para->geom->imax, jmax = para->geom->jmax;
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);  
  REAL *gx = var[GX], *gy = var[GY], *gz = var[GZ];  
  REAL *x = var[X], *y = var[Y], *z = var[Z];  
  REAL y_1,z_1;


    for(i=0;i<=imax;i++)
          for(j=0;j<=jmax;j++)
       for(k=0;k<=kmax;k++)
        {
          if(i==0 || i==imax) v[FIX(i,j,k)]= (v0[FIX(i,j,k)]+v0[FIX(i,j,k+1)])/2.0f;
        else if (k==0 || k==kmax) v[FIX(i,j,k)]= (v0[FIX(i,j,k)]+v0[FIX(i+1,j,k)])/2.0f;
        else 
        {
           y_1= (gx[FIX(i,j,k)]-x[FIX(i,j,k)])/(x[FIX(i+1,j,k)]-x[FIX(i,j,k)]);
                   z_1= (gz[FIX(i,j,k)]-z[FIX(i,j,k)])/(z[FIX(i,j,k+1)]-z[FIX(i,j,k)]);
           v[FIX(i,j,k)]= interpolation_linear(y_1, z_1, v0[FIX(i,j,k)], v0[FIX(i,j,k+1)],
                                                         v0[FIX(i+1,j,k)], v0[FIX(i+1,j,k+1)]);
        
        }
        }
}// End of v_vortex()

///////////////////////////////////////////////////////////////////////////////
///\brief  Interpolate W from surface cell to grid point
///
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param w   Pointer to W at grid point
///\param w0   Pointer to W at surface cell
///
///\return void No return needed
///////////////////////////////////////////////////////////////////////////////
void w_vortex(PARA_DATA *para, REAL **var, REAL *w, REAL *w0)
{
  int i, j, k;
  int imax = para->geom->imax, jmax = para->geom->jmax;
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);  
  REAL *gx = var[GX], *gy = var[GY], *gz = var[GZ];  
  REAL *x = var[X], *y = var[Y], *z = var[Z];  
  REAL y_1,z_1;


    for(i=0;i<=imax;i++)
          for(j=0;j<=jmax;j++)
       for(k=0;k<=kmax;k++)
        {
          if(j==0 || j==jmax) w[FIX(i,j,k)]= (w0[FIX(i,j,k)]+w0[FIX(i+1,j,k)])/2.0f;
        else if (i==0 || i==imax) w[FIX(i,j,k)]= (w0[FIX(i,j,k)]+w0[FIX(i,j+1,k)])/2.0f;
        else 
        {
           y_1= (gy[FIX(i,j,k)]-y[FIX(i,j,k)])/(y[FIX(i,j+1,k)]-y[FIX(i,j,k)]);
                   z_1= (gx[FIX(i,j,k)]-x[FIX(i,j,k)])/(x[FIX(i+1,j,k)]-x[FIX(i,j,k)]);
           w[FIX(i,j,k)]= interpolation_linear(y_1, z_1, w0[FIX(i,j,k)], w0[FIX(i+1,j,k)],
                                                         w0[FIX(i,j+1,k)], w0[FIX(i+1,j+1,k)]);
        
        }
        }

}// End of w_vortex()

///////////////////////////////////////////////////////////////////////////////
///\brief  Interpolate scalors from cell center to grid point
///
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param p   Pointer to scalor at grid point
///\param p0   Pointer to scalor at cell center
///
///\return void No return needed
///////////////////////////////////////////////////////////////////////////////
void p_vortex(PARA_DATA *para, REAL **var, REAL *p, REAL *p0)
{
  int i, j, k;
  int imax = para->geom->imax, jmax = para->geom->jmax;
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);  
  REAL *gx = var[GX], *gy = var[GY], *gz = var[GZ];  
  REAL *x = var[X], *y = var[Y], *z = var[Z];  
  REAL x_1,y_1,z_1;


    for(i=1;i<imax;i++)
          for(j=1;j<jmax;j++)
       for(k=1;k<kmax;k++)
        {
                   x_1= (gx[FIX(i,j,k)]-x[FIX(i,j,k)])/(x[FIX(i+1,j,k)]-x[FIX(i,j,k)]);
           y_1= (gy[FIX(i,j,k)]-y[FIX(i,j,k)])/(y[FIX(i,j+1,k)]-y[FIX(i,j,k)]);
                   z_1= (gz[FIX(i,j,k)]-z[FIX(i,j,k)])/(z[FIX(i,j,k+1)]-z[FIX(i,j,k)]);
           p[FIX(i,j,k)]= interpolation_bilinear(x_1, y_1, z_1, 
                                             p0[FIX(i,j,k)],  p0[FIX(i,j+1,k)], 
                                             p0[FIX(i+1,j,k)],  p0[FIX(i+1,j+1,k)], 
                                             p0[FIX(i,j,k+1)],p0[FIX(i,j+1,k+1)],
                                             p0[FIX(i+1,j,k+1)],p0[FIX(i+1,j+1,k+1)]);
        
        }
        

    i=0;
    for(j=1;j<jmax;j++)
       for(k=1;k<kmax;k++)
       {
           y_1= (gy[FIX(i,j,k)]-y[FIX(i,j,k)])/(y[FIX(i,j+1,k)]-y[FIX(i,j,k)]);
                   z_1= (gz[FIX(i,j,k)]-z[FIX(i,j,k)])/(z[FIX(i,j,k+1)]-z[FIX(i,j,k)]);
           p[FIX(i,j,k)]= interpolation_linear(y_1, z_1, p0[FIX(i,j,k)], p0[FIX(i,j,k+1)],
                                                         p0[FIX(i,j+1,k)], p0[FIX(i,j+1,k+1)]);
              
         }

    i=imax;
    for(j=1;j<jmax;j++)
       for(k=1;k<kmax;k++)
       {
           y_1= (gy[FIX(i,j,k)]-y[FIX(i,j,k)])/(y[FIX(i,j+1,k)]-y[FIX(i,j,k)]);
                   z_1= (gz[FIX(i,j,k)]-z[FIX(i,j,k)])/(z[FIX(i,j,k+1)]-z[FIX(i,j,k)]);
           p[FIX(i,j,k)]= interpolation_linear(y_1, z_1, p0[FIX(i,j,k)], p0[FIX(i,j,k+1)],
                                                         p0[FIX(i,j+1,k)], p0[FIX(i,j+1,k+1)]);
              
         }

  j=0;
      for(i=1;i<imax;i++)
       for(k=1;k<kmax;k++)
       {
             y_1= (gx[FIX(i,j,k)]-x[FIX(i,j,k)])/(x[FIX(i+1,j,k)]-x[FIX(i,j,k)]);
                   z_1= (gz[FIX(i,j,k)]-z[FIX(i,j,k)])/(z[FIX(i,j,k+1)]-z[FIX(i,j,k)]);
           p[FIX(i,j,k)]= interpolation_linear(y_1, z_1, p0[FIX(i,j,k)], p0[FIX(i,j,k+1)], 
                                                         p0[FIX(i+1,j,k)], p0[FIX(i+1,j,k+1)]);
       }
  j=jmax;
      for(i=1;i<imax;i++)
       for(k=1;k<kmax;k++)
       {
             y_1= (gx[FIX(i,j,k)]-x[FIX(i,j,k)])/(x[FIX(i+1,j,k)]-x[FIX(i,j,k)]);
                   z_1= (gz[FIX(i,j,k)]-z[FIX(i,j,k)])/(z[FIX(i,j,k+1)]-z[FIX(i,j,k)]);
           p[FIX(i,j,k)]= interpolation_linear(y_1, z_1, p0[FIX(i,j,k)], p0[FIX(i,j,k+1)],
                                                         p0[FIX(i+1,j,k)], p0[FIX(i+1,j,k+1)]);
       }
     k=0;
       for(i=1;i<imax;i++)
          for(j=1;j<jmax;j++)
        {
           y_1= (gy[FIX(i,j,k)]-y[FIX(i,j,k)])/(y[FIX(i,j+1,k)]-y[FIX(i,j,k)]);
                   z_1= (gx[FIX(i,j,k)]-x[FIX(i,j,k)])/(x[FIX(i+1,j,k)]-x[FIX(i,j,k)]);
           p[FIX(i,j,k)]= interpolation_linear(y_1, z_1, p0[FIX(i,j,k)], p0[FIX(i+1,j,k)], 
                                                         p0[FIX(i,j+1,k)], p0[FIX(i+1,j+1,k)]);
        
        }
     k=kmax;
       for(i=1;i<imax;i++)
          for(j=1;j<jmax;j++)
        {
           y_1= (gy[FIX(i,j,k)]-y[FIX(i,j,k)])/(y[FIX(i,j+1,k)]-y[FIX(i,j,k)]);
                   z_1= (gx[FIX(i,j,k)]-x[FIX(i,j,k)])/(x[FIX(i+1,j,k)]-x[FIX(i,j,k)]);
           p[FIX(i,j,k)]= interpolation_linear(y_1, z_1, p0[FIX(i,j,k)], p0[FIX(i+1,j,k)],
                                                         p0[FIX(i,j+1,k)], p0[FIX(i+1,j+1,k)]);
        
        }

}// End of p_vortex()


///////////////////////////////////////////////////////////////////////////////
///\brief  Calculate the corner value of variables
///
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param psi   Pointer to the variable
///
///\return void No return needed
///////////////////////////////////////////////////////////////////////////////
void corners(PARA_DATA *para, REAL **var, REAL *psi) {
  int i,j,k;
  int imax=para->geom->imax, jmax=para->geom->jmax;
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);



  for(k=1;k<=kmax;k++){
    psi[FIX(0,0,k)]   = (psi[FIX(1,0,k)]+psi[FIX(0,1,k)])/2.0f;
    psi[FIX(0,jmax+1,k)]= (psi[FIX(1,jmax+1,k)]+psi[FIX(0,jmax,k)])/2.0f;
    psi[FIX(imax+1,0,k)]= (psi[FIX(imax,0,k)]+psi[FIX(imax+1,1,k)])/2.0f;
    psi[FIX(imax+1,jmax+1,k)]=(psi[FIX(imax,jmax+1,k)]+psi[FIX(imax+1,jmax,k)])/2.0f;
  }

  for(j=1;j<=jmax;j++){
    psi[FIX(0,j,0)]   = (psi[FIX(1,j,0)]+psi[FIX(0,j,1)])/2.0f;
    psi[FIX(0,j,kmax+1)]= (psi[FIX(1,j,kmax+1)]+psi[FIX(0,j,kmax)])/2.0f;
    psi[FIX(imax+1,j,0)]= (psi[FIX(imax,j,0)]+psi[FIX(imax+1,j,1)])/2.0f;
    psi[FIX(imax+1,j,kmax+1)]=(psi[FIX(imax,j,kmax+1)]+psi[FIX(imax+1,j,kmax)])/2.0f;
  }

  for(i=1;i<=imax;i++){
    psi[FIX(i,0,0)]   = (psi[FIX(i,1,0)]+psi[FIX(i,0,1)])/2.0f;
    psi[FIX(i,0,kmax+1)]= (psi[FIX(i,1,kmax+1)]+psi[FIX(i,0,kmax)])/2.0f;
    psi[FIX(i,jmax+1,0)]= (psi[FIX(i,jmax,0)]+psi[FIX(i,jmax+1,1)])/2.0f;
    psi[FIX(i,jmax+1,kmax+1)]=(psi[FIX(i,jmax,kmax+1)]+psi[FIX(i,jmax+1,kmax)])/2.0f;
  }

  //W-S-B
  psi[FIX(0,0,0)] = (psi[FIX(0,1,0)]+psi[FIX(1,0,0)]+psi[FIX(0,0,1)]) / 3.0f;

  //W-N-B
  psi[FIX(0,jmax+1,0)] = ( psi[FIX(1,jmax+1,0)]+psi[FIX(0,jmax,0)]
                         +psi[FIX(0,jmax+1,1)]) / 3.0f;
  //E-S-B
  psi[FIX(imax+1,0,0)] = ( psi[FIX(imax,0,0)]+psi[FIX(imax+1,1,0)]
                         +psi[FIX(imax+1,0,1)]) / 3.0f;
  //E-N-B
  psi[FIX(imax+1,jmax+1,0)] = ( psi[FIX(imax,jmax+1,0)]+psi[FIX(imax+1,jmax,0)]
                              +psi[FIX(imax+1,jmax+1,1)]) / 3.0f;
  //W-S-F
  psi[FIX(0,0,kmax+1)] = ( psi[FIX(0,1,kmax+1)]+psi[FIX(1,0,kmax+1)]
                         +psi[FIX(0,0,kmax)]) / 3.0f;  
  //W-N-F
  psi[FIX(0,jmax+1,kmax+1)] = ( psi[FIX(1,jmax+1,kmax+1)]+psi[FIX(0,jmax,kmax+1)]
                              +psi[FIX(0,jmax+1,kmax)]) / 3.0f;

  //E-S-F
  psi[FIX(imax+1,0,kmax+1)] = ( psi[FIX(imax,0,kmax+1)]+psi[FIX(imax+1,1,kmax+1)]
                              +psi[FIX(imax+1,0,kmax)]) / 3.0f;
  //E-N-F
  psi[FIX(imax+1,jmax+1,kmax+1)] = ( psi[FIX(imax,jmax+1,kmax+1)]+psi[FIX(imax+1,jmax,kmax+1)]
                                   +psi[FIX(imax+1,jmax+1,kmax)]) / 3.0f;


} // End of corners()


///////////////////////////////////////////////////////////////////////////////
///\brief  Calculate the edge value of variables
///
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param psi   Pointer to the variable
///
///\return void No return needed
///////////////////////////////////////////////////////////////////////////////
void corners_vortex(PARA_DATA *para, REAL **var, REAL *psi)
{
  int i,j,k;
  int imax=para->geom->imax, jmax=para->geom->jmax;
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);



   for(k=1;k<kmax;k++)
    {
		psi[FIX(0,0,k)]   = (psi[FIX(1,0,k)]+psi[FIX(0,1,k)])/2.0f;
		psi[FIX(0,jmax,k)]= (psi[FIX(1,jmax,k)]+psi[FIX(0,jmax-1,k)])/2.0f;
		psi[FIX(imax,0,k)]= (psi[FIX(imax-1,0,k)]+psi[FIX(imax,1,k)])/2.0f;
		psi[FIX(imax,jmax,k)]=(psi[FIX(imax-1,jmax,k)]+psi[FIX(imax,jmax-1,k)])/2.0f;
    }

  for(j=1;j<jmax;j++)
    {
		psi[FIX(0,j,0)]   = (psi[FIX(1,j,0)]+psi[FIX(0,j,1)])/2.0f;
		psi[FIX(0,j,kmax)]= (psi[FIX(1,j,kmax)]+psi[FIX(0,j,kmax-1)])/2.0f;
		psi[FIX(imax,j,0)]= (psi[FIX(imax-1,j,0)]+psi[FIX(imax,j,1)])/2.0f;
		psi[FIX(imax,j,kmax)]=(psi[FIX(imax-1,j,kmax)]+psi[FIX(imax,j,kmax-1)])/2.0f;
    }

    for(i=1;i<imax;i++)
    {
		psi[FIX(i,0,0)]   = (psi[FIX(i,1,0)]+psi[FIX(i,0,1)])/2.0f;
		psi[FIX(i,0,kmax)]= (psi[FIX(i,1,kmax)]+psi[FIX(i,0,kmax-1)])/2.0f;
		psi[FIX(i,jmax,0)]= (psi[FIX(i,jmax-1,0)]+psi[FIX(i,jmax,1)])/2.0f;
		psi[FIX(i,jmax,kmax)]=(psi[FIX(i,jmax-1,kmax)]+psi[FIX(i,jmax,kmax-1)])/2.0f;
    }

  //W-S-B
  psi[FIX(0,0,0)] = (psi[FIX(0,1,0)]+psi[FIX(1,0,0)]+psi[FIX(0,0,1)]) / 3.0f;

  //W-N-B
  psi[FIX(0,jmax,0)] = ( psi[FIX(1,jmax,0)]+psi[FIX(0,jmax-1,0)]
                         +psi[FIX(0,jmax,1)]) / 3.0f;
  //E-S-B
  psi[FIX(imax,0,0)] = ( psi[FIX(imax-1,0,0)]+psi[FIX(imax,1,0)]
                         +psi[FIX(imax,0,1)]) / 3.0f;
  //E-N-B
  psi[FIX(imax,jmax,0)] = ( psi[FIX(imax-1,jmax,0)]+psi[FIX(imax,jmax-1,0)]
                              +psi[FIX(imax,jmax,1)]) / 3.0f;
  //W-S-F
  psi[FIX(0,0,kmax)] = ( psi[FIX(0,1,kmax)]+psi[FIX(1,0,kmax)]
                         +psi[FIX(0,0,kmax-1)]) / 3.0f;  
  //W-N-F
  psi[FIX(0,jmax,kmax)] = ( psi[FIX(1,jmax,kmax)]+psi[FIX(0,jmax-1,kmax)]
                              +psi[FIX(0,jmax,kmax-1)]) / 3.0f;

  //E-S-F
  psi[FIX(imax,0,kmax)] = ( psi[FIX(imax-1,0,kmax)]+psi[FIX(imax,1,kmax)]
                              +psi[FIX(imax,0,kmax-1)]) / 3.0f;
  //E-N-F
  psi[FIX(imax,jmax,kmax)] = ( psi[FIX(imax-1,jmax,kmax)]+psi[FIX(imax,jmax-1,kmax)]
                                   +psi[FIX(imax,jmax,kmax-1)]) / 3.0f;


} // End of corners_vortex()

///////////////////////////////////////////////////////////////////////////////
///\brief  Interpolate velocities to the cell center
///
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///
///\return void No return needed
///////////////////////////////////////////////////////////////////////////////
void uvw_center(PARA_DATA *para, REAL **var) {
  int i, j, k;
  int imax=para->geom->imax, jmax=para->geom->jmax;
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  REAL *u0 = var[VX], *v0 = var[VY], *w0 = var[VZ];
  REAL *flagp=var[FLAGP];

  for(j=0; j<=jmax+1; j++)
    for(k=0; k<=kmax+1; k++) {
      u0[FIX(imax+1,j,k)] = u0[FIX(imax,j,k)];
      for(i=imax; i>=1; i--) {
        if(flagp[FIX(i,j,k)]==1) u0[FIX(i,j,k)]=0;
        else u0[FIX(i,j,k)] = 0.5f * (u0[FIX(i,j,k)]+u0[FIX(i-1,j,k)]);
      }
    }
    
  for(i=0; i<=imax+1; i++)
    for(k=0; k<=kmax+1; k++) {
      v0[FIX(i,jmax+1,k)] = v0[FIX(i,jmax,k)];
      for(j=jmax; j>=1; j--) {
        if(flagp[FIX(i,j,k)]==1) v0[FIX(i,j,k)]=0;
        else v0[FIX(i,j,k)] = 0.5f * (v0[FIX(i,j,k)]+v0[FIX(i,j-1,k)]);
      }
    }
    
  for(i=0; i<=imax+1; i++)
    for(j=0; j<=jmax+1; j++) {
      w0[FIX(i,j,kmax+1)] = w0[FIX(i,j,kmax)];
      for(k=kmax; k>=1; k--) {
        if(flagp[FIX(i,j,k)]==1) w0[FIX(i,j,k)]=0;
        else w0[FIX(i,j,k)] = 0.5f * (w0[FIX(i,j,k)]+w0[FIX(i,j,k-1)]);
      }
    }
}// End of uvw_center()
