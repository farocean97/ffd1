///////////////////////////////////////////////////////////////////////////////
///
/// \file   interpolation.c
///
/// \brief  Subroutines for performing interpolation step for semi-Lagrangian
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
/// This file provides functions for performing interpolations for
/// semi-Lagrangian scheme.linear interpolation scheme \c interpolaiton_bilinear()
/// will be used in FFD, but if a heat source cell boudnary condition was assigned
/// in FFD,the upwind scheme will be used for this cell when plume model is not 
/// applied for the heat source cell.
///
///////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include "data_structure.h"
#include "interpolation.h"
#include "utility.h"


///////////////////////////////////////////////////////////////////////////////
/// \brief Entrance of interpolation
///
///\param para Pointer to FFD parameters
///\param d0 Pointer to the variable for interpolation
///\param x_1 Reciprocal of X-length
///\param y_1 Reciprocal of Y-length
///\param z_1 Reciprocal of Z-length
///\param p I-index of the control volume
///\param q J-index of the control volume
///\param r K-index of the control volume
///
///\return Interpolated value
///////////////////////////////////////////////////////////////////////////////
REAL interpolation(PARA_DATA *para, REAL *d0, REAL x_1, REAL y_1, REAL z_1,
                   int p, int q, int r) {
  int imax = para->geom->imax, jmax = para->geom->jmax;
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);

  return interpolation_bilinear(x_1, y_1, z_1, 
        d0[FIX(p,q,r)],  d0[FIX(p,q+1,r)],  d0[FIX(p+1,q,r)],  d0[FIX(p+1,q+1,r)], 
        d0[FIX(p,q,r+1)],d0[FIX(p,q+1,r+1)],d0[FIX(p+1,q,r+1)],d0[FIX(p+1,q+1,r+1)]);

} // End of interpolation()


///////////////////////////////////////////////////////////////////////////////
///\brief Bilinear interpolation
///
///\param x_1 Reciprocal of X-length
///\param y_1 Reciprocal of Y-length
///\param z_1 Reciprocal of Z-length
///\param d000 parameter for interpolation
///\param d010 parameter for interpolation
///\param d100 parameter for interpolation
///\param d110 parameter for interpolation
///\param d001 parameter for interpolation
///\param d011 parameter for interpolation
///\param d101 parameter for interpolation
///\param d111 parameter for interpolation
//
///\return Interpolated value
///////////////////////////////////////////////////////////////////////////////
REAL interpolation_bilinear(REAL x_1, REAL y_1, REAL z_1,
                            REAL d000, REAL d010, REAL d100, REAL d110,
                            REAL d001, REAL d011, REAL d101, REAL d111){
  REAL x_0, y_0, z_0;
  REAL tmp0, tmp1;

  x_0 = 1.0f - x_1;
  y_0 = 1.0f - y_1; 
  z_0 = 1.0f - z_1;
  
  tmp0 = x_0*(y_0*d000+y_1*d010) + x_1*(y_0*d100+y_1*d110);
  tmp1 = x_0*(y_0*d001+y_1*d011) + x_1*(y_0*d101+y_1*d111);


  return z_0*tmp0+z_1*tmp1;   

} // End of interpolation_bilinear()



///////////////////////////////////////////////////////////////////////////////
///\brief Bilinear interpolation for surface value
///
///\param x_1 Reciprocal of X-length
///\param y_1 Reciprocal of Y-length
///\param d000 parameter for interpolation
///\param d010 parameter for interpolation
///\param d100 parameter for interpolation
///\param d110 parameter for interpolation
///
///\return Interpolated value
///////////////////////////////////////////////////////////////////////////////
REAL interpolation_linear(REAL x_1, REAL y_1, REAL d000, 
                          REAL d010, REAL d100, REAL d110){
  REAL x_0, y_0;
  REAL tmp0;

  x_0 = 1.0f - x_1; 
  y_0 = 1.0f - y_1;
  
  tmp0 = x_0*(y_0*d000+y_1*d010) + x_1*(y_0*d100+y_1*d110);

  return tmp0;

} // End of interpolation_bilinear()


///////////////////////////////////////////////////////////////////////////////
///\brief Upwind advection scheme of heat source cells
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD variables
///\param d0 Pointer to the variable for interpolation
///\param i I-index of the control volume
///\param j J-index of the control volume
///\param k K-index of the control volume
///\param d110 parameter for interpolation
///
///\return advection value
///////////////////////////////////////////////////////////////////////////////
REAL advection_upwind(PARA_DATA *para, REAL **var,REAL *d0,int i, int j, int k) { 
  int imax = para->geom->imax, jmax = para->geom->jmax;
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  REAL Dx,Dy,Dz;
  REAL Fw,Fe,Fn,Fs,Ff,Fb;
  REAL *x = var[X], *y = var[Y], *z = var[Z];
  REAL *u=var[VX],*v=var[VY],*w=var[VZ];
  REAL *gx = var[GX], *gy = var[GY], *gz = var[GZ];
  REAL dt = para->mytime->dt; 
  REAL d;

  Dx=gx[FIX(i,j,k)]-gx[FIX(i-1,j,k)];
  Dy=gy[FIX(i,j,k)]-gy[FIX(i,j-1,k)];
  Dz=gz[FIX(i,j,k)]-gz[FIX(i,j,k-1)];
  Fw=Dy*Dz*u[FIX(i-1,j,k)];
  Fe=Dy*Dz*u[FIX(i,j,k)];
  Fn=Dx*Dz*v[FIX(i,j,k)];
  Fs=Dx*Dz*v[FIX(i,j-1,k)];
  Ff=Dy*Dx*w[FIX(i,j,k)];
  Fb=Dy*Dx*w[FIX(i,j,k-1)];

  
  d=d0[FIX(i,j,k)]-dt/(Dx*Dy*Dz)*( max(Fe,0)*d0[FIX(i,j,k)]-max(Fw,0)*d0[FIX(i-1,j,k)]
                               + max(Fn,0)*d0[FIX(i,j,k)]-max(Fs,0)*d0[FIX(i,j-1,k)]
                 + max(Ff,0)*d0[FIX(i,j,k)]-max(Fb,0)*d0[FIX(i,j,k-1)]
                 + min(Fe,0)*d0[FIX(i+1,j,k)]-min(Fw,0)*d0[FIX(i,j,k)]
                               + min(Fn,0)*d0[FIX(i,j+1,k)]-min(Fs,0)*d0[FIX(i,j,k)]
                 + min(Ff,0)*d0[FIX(i,j,k+1)]-min(Fb,0)*d0[FIX(i,j,k)]);

  return d;

} //End of advection_upwind()


///////////////////////////////////////////////////////////////////////////////
///\brief Modifying the interpolation coefficients of boundary cells
///
/// For near boundary cells, it's possible that the interpolation would
/// interoplate from boundary cells, which will generate unphysical 
/// results. So it is necessary to adjust the coefficients of boundary cells.
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD variables
///\param coef Pointer to interpolation coefficients of surounding cells
///\param x1 Reciprocal of X-length
///\param y1 Reciprocal of Y-length
///\param z1 Reciprocal of Z-length
///\param p I-index of the control volume
///\param q J-index of the control volume
///\param r K-index of the control volume
///
///\return advection value
/////////////////////////////////////////////////////////////////////////////// 
void interpolation_coef(PARA_DATA *para,REAL **var, REAL *coef,REAL x1, REAL y1, REAL z1,int p, int q, int r){

  int imax = para->geom->imax, jmax = para->geom->jmax;
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  REAL x0, y0, z0;
  int i,j,k;
  int ii,jj,kk;
  REAL sum=1;
  int cell[8]={0};
  REAL *flagp = var[FLAGP],*flagu = var[FLAGU];
  REAL *flagv = var[FLAGV],*flagw = var[FLAGW];
  REAL *flagt=var[FLAGT];
  REAL eta;
  REAL total;

  x0 = 1.0f - x1;
  y0 = 1.0f - y1; 
  z0 = 1.0f - z1;
 
  coef[PIX(0,0,0)]=x0*y0*z0;
  coef[PIX(0,1,0)]=x0*y1*z0;
  coef[PIX(1,0,0)]=x1*y0*z0;
  coef[PIX(1,1,0)]=x1*y1*z0;

  coef[PIX(0,0,1)]=x0*y0*z1;
  coef[PIX(0,1,1)]=x0*y1*z1;
  coef[PIX(1,0,1)]=x1*y0*z1;
  coef[PIX(1,1,1)]=x1*y1*z1;

 
  for(i=0;i<=1;i++)
    for(j=0;j<=1;j++)
      for(k=0;k<=1;k++) {
        if( flagp[FIX(p+i,q+j,r+k)]>0  && flagt[FIX(p+i,q+j,r+k)]<0 && coef[PIX(i,j,k)]<1.0) {
          eta=sum/(sum-coef[PIX(i,j,k)]);
          coef[PIX(i,j,k)]=0;
          total=0;
          for(ii=0;ii<=1;ii++)
            for(jj=0;jj<=1;jj++)
              for(kk=0;kk<=1;kk++){
                coef[PIX(ii,jj,kk)] *= eta;
                total += coef[PIX(ii,jj,kk)];
              }
        }
      }
 } // End of interpolation_bilinear()



///////////////////////////////////////////////////////////////////////////////
///\brief Preforming interpolation for temperature and specices.
///
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD variables
///\param d0 Pointer to interpolation variable
///\param x_1 Reciprocal of X-length
///\param y_1 Reciprocal of Y-length
///\param z_1 Reciprocal of Z-length
///\param i0 I-index of the arriving point
///\param j0 J-index of the arriving point
///\param k0 K-index of the carriving point
///\param p I-index of the control volume
///\param q J-index of the control volume
///\param r K-index of the control volume
///
///\return advection value
/////////////////////////////////////////////////////////////////////////////// 
REAL interpolation_temp(PARA_DATA *para, REAL **var,REAL *d0, REAL x_1, REAL y_1, REAL z_1,
                   int i0, int j0, int k0, int p, int q, int r) {
  int imax = para->geom->imax, jmax = para->geom->jmax;
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  int i,j,k;
  int ii,jj,kk;
  REAL coef[8];
  REAL d[8];
  REAL interp_value=0;
  REAL tcoef=0;
  REAL *flagu = var[FLAGU];
  REAL *flagv = var[FLAGV],*flagw = var[FLAGW];
  REAL *x = var[X],  *y = var[Y],  *z = var[Z]; 
  REAL *gx = var[GX],  *gy = var[GY],  *gz = var[GZ]; 
  REAL lx,ly,lz;
  
  interpolation_coef(para,var, coef,x_1, y_1, z_1,p,q,r);

  for(i=0;i<=1;i++)
    for(j=0;j<=1;j++)
      for(k=0;k<=1;k++){
        d[PIX(i,j,k)]=d0[FIX(p+i,q+j,r+k)];
      }

  lx=x_1*(x[FIX(p+1,q,r)]-x[FIX(p,q,r)]);
  ly=y_1*(y[FIX(p,q+1,r)]-y[FIX(p,q,r)]);
  lz=z_1*(z[FIX(p,q,r+1)]-z[FIX(p,q,r)]);

  if(lx<=(gx[FIX(p,q,r)]-x[FIX(p,q,r)])) ii=1;
  else ii=0;

  if(ly<=(gy[FIX(p,q,r)]-y[FIX(p,q,r)])) jj=1;
  else jj=0;

  if(lz<=(gz[FIX(p,q,r)]-z[FIX(p,q,r)])) kk=1;
  else kk=0;
  
  //Special treatment for cells around a partition wall
  if(flagu[FIX(p,q+1-jj,r+1-kk)]==4)     d[PIX(ii,1-jj,1-kk)]=d0[FIX(p+1-ii,q+1-jj,r+1-kk)];
  if(flagv[FIX(p+1-ii,q,r+1-kk)]==4)     d[PIX(1-ii,jj,1-kk)]=d0[FIX(p+1-ii,q+1-jj,r+1-kk)];
  if(flagu[FIX(p,q+jj,r+1-kk)]==4)       d[PIX(ii,jj,1-kk)]=d0[FIX(p+1-ii,q+jj,r+1-kk)];
  if(flagv[FIX(p+ii,q,r+1-kk)]==4)       d[PIX(ii,jj,1-kk)]=d0[FIX(p+ii,q+1-jj,r+1-kk)];
  if(flagw[FIX(p+1-ii,q+1-jj,r)]==4)     d[PIX(1-ii,1-jj,kk)]=d0[FIX(p+1-ii,q+1-jj,r+1-kk)];
  if(flagw[FIX(p+1-ii,q+jj,r)]==4)       d[PIX(1-ii,jj,kk)]=d0[FIX(p+1-ii,q+jj,r+1-kk)];
  if(flagw[FIX(p+ii,q+1-jj,r)]==4)       d[PIX(ii,1-jj,kk)]=d0[FIX(p+ii,q+1-jj,r+1-kk)];
  if(flagw[FIX(p+ii,q+jj,r)]==4)         d[PIX(ii,jj,kk)]=d0[FIX(p+ii,q+jj,r+1-kk)];
  if(flagu[FIX(p,q+1-jj,r+kk)]==4)       d[PIX(ii,1-jj,kk)]=d0[FIX(p+1-ii,q+1-jj,r+kk)];
  if(flagv[FIX(p+1-ii,q,r+kk)]==4)       d[PIX(1-ii,jj,kk)]=d0[FIX(p+1-ii,q+1-jj,r+kk)];
  if(flagu[FIX(p,q+jj,r+kk)]==4)         d[PIX(ii,jj,kk)]=d0[FIX(p+1-ii,q+jj,r+kk)];
  if(flagv[FIX(p+ii,q,r+kk)]==4)         d[PIX(ii,jj,kk)]=d0[FIX(p+ii,q+1-jj,r+kk)];

  
  for(i=0;i<=1;i++)
    for(j=0;j<=1;j++)
      for(k=0;k<=1;k++){
        interp_value += coef[PIX(i,j,k)]*d[PIX(i,j,k)];
        tcoef +=coef[PIX(i,j,k)];
      }

   var[LOCMIN][FIX(i0,j0,k0)]=check_min_pix(para, d); 
   var[LOCMAX][FIX(i0,j0,k0)]=check_max_pix(para, d);

  return interp_value;

} // End of interpolation_temp()

