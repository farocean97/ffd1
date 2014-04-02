#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "data_structure.h"
#include "utility.h"
#include "interpolation.h"


REAL check_residual(PARA_DATA *para, REAL **var, REAL *flag,REAL *x)
{
  int imax = para->geom->imax, jmax = para->geom->jmax; 
  int kmax = para->geom->kmax;
  int i, j, k;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  REAL *aw = var[AW], *ae = var[AE], *as = var[AS], *an = var[AN];
  REAL *ap = var[AP], *ab = var[AB], *af = var[AF], *b = var[B];  
  REAL tmp1=0,tmp2=0, residual = 0.0; 



  FOR_EACH_CELL
    if(flag[FIX(i,j,k)]>=0) continue;
    tmp1 += (float) fabs(ap[FIX(i,j,k)]*x[FIX(i,j,k)] 
        - ae[FIX(i,j,k)]*x[FIX(i+1,j,k)] - aw[FIX(i,j,k)]*x[FIX(i-1,j,k)]
        - an[FIX(i,j,k)]*x[FIX(i,j+1,k)] - as[FIX(i,j,k)]*x[FIX(i,j-1,k)]
        - af[FIX(i,j,k)]*x[FIX(i,j,k+1)] - ab[FIX(i,j,k)]*x[FIX(i,j,k-1)]
        - b[FIX(i,j,k)]);
    tmp2 += (float) fabs(ap[FIX(i,j,k)]*x[FIX(i,j,k)]);
  END_FOR
    
  return tmp1/tmp2;

}// End of check_residual( )


REAL check_integ(PARA_DATA *para, REAL **var)
{
  int imax = para->geom->imax, jmax = para->geom->jmax; 
  int kmax = para->geom->kmax;
  int i, j, k;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  REAL *u=var[VX],*v=var[VY],*w=var[VZ];
  REAL *gx=var[GX],*gy=var[GY],*gz=var[GZ];
  REAL *flag=var[FLAGP];
  REAL uu,vv,ww;
  REAL dx,dy,dz;
  REAL vel=0;
  

  FOR_EACH_CELL
    uu=0.5f*(u[FIX(i,j,k)]+u[FIX(i-1,j,k)]);
      vv=0.5f*(v[FIX(i,j,k)]+v[FIX(i,j-1,k)]);
    ww=0.5f*(w[FIX(i,j,k)]+w[FIX(i,j,k-1)]);
    dx=gx[FIX(i,j,k)]-gx[FIX(i-1,j,k)];
    dy=gy[FIX(i,j,k)]-gy[FIX(i,j-1,k)];
    dz=gz[FIX(i,j,k)]-gz[FIX(i,j,k-1)];

     vel += dx*dy*dz*sqrt(uu*uu+vv*vv+ww*ww);

  END_FOR
    
  return vel;

}// End of check_residual( )


void swap(PARA_DATA *para, REAL **var)
{
  int i, size;
  size=(para->geom->imax+2) * (para->geom->jmax+2) * (para->geom->kmax+2);
  
  
   for(i=0; i<size; i++) 
  {
    var[VX][i]     = var[TMP1][i];
     var[VY][i]     = var[TMP2][i];
     var[VZ][i]     = var[TMP3][i];
   //  var[TMP1][i]=var[VX][i]    ;
  //  var[TMP2][i]= var[VY][i] ;
   //  var[TMP3][i]=var[VZ][i];
  }
}


void swapUVW_C(PARA_DATA *para, REAL **var)
{
  int i, j, k;
  int imax=para->geom->imax, jmax=para->geom->jmax;
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  REAL *u = var[VX], *v = var[VY], *w = var[VZ];
  REAL *us = var[VXC], *vs = var[VYC], *ws = var[VZC];

    for(j=0; j<=jmax+1; j++)
    {
    for(k=0; k<=kmax+1; k++)
          {
            
                us[FIX(imax+1,j,k)] = u[FIX(imax,j,k)];
                us[FIX(0,j,k)] = u[FIX(0,j,k)];
                     for(i=imax; i>=1; i--)
                     {
                          us[FIX(i,j,k)] = 0.5f * (u[FIX(i,j,k)]+u[FIX(i-1,j,k)]);
                      }
            }
    }


  
    for(i=0; i<=imax+1; i++)
    {
    for(k=0; k<=kmax+1; k++)
          {
            
                vs[FIX(i,jmax+1,k)] = v[FIX(i,jmax,k)];
                vs[FIX(i,0,k)] = v[FIX(i,0,k)];
                     for(j=jmax; j>=1; j--)
                     {
                          vs[FIX(i,j,k)] = 0.5f * (v[FIX(i,j,k)]+v[FIX(i,j-1,k)]);
                      }
            }
    }
   

   for(i=0; i<=imax+1; i++)
    {
    for(j=0; j<=jmax+1; j++)
          {
            
                ws[FIX(i,j,kmax+1)] = w[FIX(i,j,kmax)];
                ws[FIX(i,j,0)] = w[FIX(i,j,0)];

                     for(k=kmax; k>=1; k--)
                     {
                          ws[FIX(i,j,k)] = 0.5f * (w[FIX(i,j,k)]+w[FIX(i,j,k-1)]);
                      }
            }
    }


}


void swapUVW(PARA_DATA *para, REAL **var)
{
  int i, j, k;
  int imax=para->geom->imax, jmax=para->geom->jmax;
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  REAL fac;
  REAL *u = var[TMP1], *v = var[TMP2], *w = var[TMP3];
  REAL *us = var[VXS], *vs = var[VYS], *ws = var[VZS];
  REAL *x = var[X], *y = var[Y], *z =var[Z];
  REAL *gx = var[GX], *gy = var[GY], *gz =var[GZ];

for(i=1; i<=imax-1; i++)
     for(j=1; j<=jmax; j++)
    for(k=1; k<=kmax; k++)
          {
         fac=(gx[FIX(i,j,k)]-x[FIX(i,j,k)])/(x[FIX(i+1,j,k)]-x[FIX(i,j,k)]);
               u[FIX(i,j,k)]=(1-fac)*us[FIX(i,j,k)]+fac*us[FIX(i+1,j,k)];
      }

for(i=1; i<=imax; i++)
     for(j=1; j<=jmax-1; j++)
    for(k=1; k<=kmax; k++)
          {
         fac=(gy[FIX(i,j,k)]-y[FIX(i,j,k)])/(y[FIX(i,j+1,k)]-y[FIX(i,j,k)]);
               v[FIX(i,j,k)]=(1-fac)*vs[FIX(i,j,k)]+fac*vs[FIX(i,j+1,k)];
      }

for(i=1; i<=imax; i++)
     for(j=1; j<=jmax; j++)
    for(k=1; k<=kmax-1; k++)
          {
         fac=(gz[FIX(i,j,k)]-z[FIX(i,j,k)])/(z[FIX(i,j,k+1)]-z[FIX(i,j,k)]);
               w[FIX(i,j,k)]=(1-fac)*ws[FIX(i,j,k)]+fac*ws[FIX(i,j,k+1)];
      }

}


void swapT(PARA_DATA *para, REAL **var)
{
  int i, j, k;
  int imax=para->geom->imax, jmax=para->geom->jmax;
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  REAL *flagp=var[FLAGP];

  for(i=1; i<=imax; i++)
     for(j=1; j<=jmax; j++)
    for(k=1; k<=kmax; k++)
    {
         if(flagp[FIX(i,j,k)]==1) continue;
             var[TEMP][FIX(i,j,k)]     = var[TMP1][FIX(i,j,k)];              
        
               } 

}

void swapT0(PARA_DATA *para, REAL **var)
{
  int i, j, k;
  int imax=para->geom->imax, jmax=para->geom->jmax;
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  REAL *flagp=var[FLAGP];

  for(i=1; i<=imax; i++)
     for(j=1; j<=jmax; j++)
    for(k=1; k<=kmax; k++)
    {
         if(flagp[FIX(i,j,k)]==1) continue;
             var[TMP1][FIX(i,j,k)]     = var[TEMP][FIX(i,j,k)];              
        
               } 

}


void limit(REAL x)
{

  REAL min,max;

  min=-1e6;
  max=1e6;
  x=x>min?x:min;
  x=x<max?x:max;

}

void check_mass(PARA_DATA *para, REAL **var,int i, int j,int k)
{
  int imax = para->geom->imax, jmax = para->geom->jmax; 
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  REAL *u=var[VX],*v=var[VY],*w=var[VZ],*gx=var[GX],*gy=var[GY],*gz=var[GZ];
  REAL mass=0;
 

   mass  = (u[FIX(i,j,k)]-u[FIX(i-1,j,k)])*(gy[FIX(i,j,k)]-gy[FIX(i,j-1,k)])*(gz[FIX(i,j,k)]-gz[FIX(i,j,k-1)])
      + (v[FIX(i,j,k)]-v[FIX(i,j-1,k)])*(gx[FIX(i,j,k)]-gx[FIX(i-1,j,k)])*(gz[FIX(i,j,k)]-gz[FIX(i,j,k-1)])
      + (w[FIX(i,j,k)]-w[FIX(i,j,k-1)])*(gx[FIX(i,j,k)]-gx[FIX(i-1,j,k)])*(gy[FIX(i,j,k)]-gy[FIX(i,j-1,k)]);

  mass = mass/((v[FIX(i,j,k)]-v[FIX(i,j-1,k)])*(gx[FIX(i,j,k)]-gx[FIX(i-1,j,k)])*(gz[FIX(i,j,k)]-gz[FIX(i,j,k-1)]));

  //printf("mass(%d,%d,%d)=%f\n",i,j,k,mass);

 
}// End of check_residual( )

REAL check_mass_W(PARA_DATA *para, REAL **var,int i1, int i2,int j1,int j2,int k1,int k2)
{
  int imax = para->geom->imax, jmax = para->geom->jmax; 
  int kmax = para->geom->kmax;
  int i, j, k;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  REAL *w=var[VZ],*gx=var[GX],*gy=var[GY],*gz=var[GZ];
  REAL mass=0;
 
 for(i=i1;i<=i2;i++)
  for(j=j1;j<=j2;j++)
    for(k=k1;k<=k2;k++)
    {
      mass  += w[FIX(i,j,k)]*(gy[FIX(i,j,k)]-gy[FIX(i,j-1,k)])*(gx[FIX(i-1,j,k)]-gx[FIX(i,j,k)]);
    }

    return mass;
 
}// End of check_residual( )

REAL check_mass_V(PARA_DATA *para, REAL **var,int i1, int i2,int j1,int j2,int k1,int k2)
{
  int imax = para->geom->imax, jmax = para->geom->jmax; 
  int kmax = para->geom->kmax;
  int i, j, k;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  REAL *v=var[VY],*gx=var[GX],*gy=var[GY],*gz=var[GZ];
  REAL mass=0;
 
 for(i=i1;i<=i2;i++)
  for(j=j1;j<=j2;j++)
    for(k=k1;k<=k2;k++)
    {
      mass  += v[FIX(i,j,k)]*(gz[FIX(i,j,k)]-gz[FIX(i,j,k-1)])*(gx[FIX(i,j,k)]-gx[FIX(i-1,j,k)]);
     // printf("v=%f\n",v[FIX(i,j,k)]);
    }

    return mass;
 
}// End of check_residual( )


REAL check_mass_in(PARA_DATA *para, REAL **var,int i1, int i2,int j1,int j2,int k1,int k2)
{
  int imax = para->geom->imax, jmax = para->geom->jmax; 
  int kmax = para->geom->kmax;
  int i, j, k;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  REAL *u=var[VX],*gx=var[GX],*gy=var[GY],*gz=var[GZ];
  REAL mass=0;
 
 for(i=i1;i<=i2;i++)
  for(j=j1;j<=j2;j++)
    for(k=k1;k<=k2;k++)
    {
      if (u[FIX(i,j,k)]>0)
      {
      mass  += u[FIX(i,j,k)]*(gy[FIX(i,j,k)]-gy[FIX(i,j-1,k)])*(gz[FIX(i,j,k)]-gz[FIX(i,j,k-1)]);
      }
    }

    return mass;
 
}// End of check_residual( )


REAL check_mass_out(PARA_DATA *para, REAL **var,int i1, int i2,int j1,int j2,int k1,int k2)
{
  int imax = para->geom->imax, jmax = para->geom->jmax; 
  int kmax = para->geom->kmax;
  int i, j, k;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  REAL *u=var[VX],*gx=var[GX],*gy=var[GY],*gz=var[GZ];
  REAL mass=0;
 
 for(i=i1;i<=i2;i++)
  for(j=j1;j<=j2;j++)
    for(k=k1;k<=k2;k++)
    {
       mass  += u[FIX(i,j,k)]*(gy[FIX(i,j,k)]-gy[FIX(i,j-1,k)])*(gz[FIX(i,j,k)]-gz[FIX(i,j,k-1)]);
    }

    return mass;
 
}// End of check_residual( )


REAL check_particle(PARA_DATA *para, REAL **var,int i1, int i2,int j1,int j2,int k1,int k2)
{
  int imax = para->geom->imax, jmax = para->geom->jmax; 
  int kmax = para->geom->kmax;
  int i, j, k,it;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  REAL *u=var[VX],*gx=var[GX],*gy=var[GY],*gz=var[GZ];
  REAL *den=var[DEN];
  REAL mass=0;
 
 for(i=i1;i<=i2;i++)
  for(j=j1;j<=j2;j++)
    for(k=k1;k<=k2;k++)
    {
      mass  += den[FIX(i,j,k)]*u[FIX(i,j,k)]*(gy[FIX(i,j,k)]-gy[FIX(i,j-1,k)])*(gz[FIX(i,j,k)]-gz[FIX(i,j,k-1)]);
    }

    return mass;
 
}// End of check_residual( )


void psi_conservation(PARA_DATA *para, REAL **var, REAL *psi,REAL *psi0,int **BINDEX) {
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
        var[LOCMIN][FIX(i,j,k)]=check_min(para, psi0, BINDEX[4][FIX(i,j,k)],  BINDEX[5][FIX(i,j,k)],  BINDEX[6][FIX(i,j,k)]); 
        var[LOCMAX][FIX(i,j,k)]=check_max(para, psi0, BINDEX[4][FIX(i,j,k)],  BINDEX[5][FIX(i,j,k)],  BINDEX[6][FIX(i,j,k)]); 
      }
    END_FOR

    qin= inflow(para,var,psi0,BINDEX);
    qout=outflow(para,var,psi0,BINDEX);
    mass0 +=(qin+qout)*dt;

    FOR_EACH_CELL
      if(flagp[FIX(i,j,k)]>0) continue;
      dA= (gx[FIX(i,j,k)]-gx[FIX(i-1,j,k)])*(gz[FIX(i,j,k)]-gz[FIX(i,j,k-1)])*(gy[FIX(i,j,k)]-gy[FIX(i,j-1,k)]);
      area += dA;
      if(BINDEX[20][FIX(i,j,k)]==1) mass0 += rho*cp*var[TMP2][FIX(i,j,k)]*dA;
      else mass0 += rho*cp*psi0[FIX(i,j,k)]*dA;
      massstar += rho*cp*psi0[FIX(i,j,k)]*dA;
      mass  += rho*cp*psi[FIX(i,j,k)]*dA;
      dens +=  (float) fabs(psi[FIX(i,j,k)]-var[LOCMIN][FIX(i,j,k)]);
      dens1 +=  (float) fabs(psi[FIX(i,j,k)]-var[LOCMAX][FIX(i,j,k)]);   
    END_FOR

    eta=mass0-mass;
    //printf("eta=%f\n",eta);
    if(eta<0) {   // mass generation occured
      FOR_EACH_CELL
        if(flagp[FIX(i,j,k)]>0 ) continue;
        dens_s[FIX(i,j,k)]=(psi[FIX(i,j,k)]-var[LOCMIN][FIX(i,j,k)])/dens*eta/(rho*cp*(gx[FIX(i,j,k)]-gx[FIX(i-1,j,k)])*(gz[FIX(i,j,k)]-gz[FIX(i,j,k-1)])*(gy[FIX(i,j,k)]-gy[FIX(i,j-1,k)]));
        //if(i==9 && j==5 && k==8) printf("dens_s=%f\n",dens_s[FIX(i,j,k)]);
        psi[FIX(i,j,k)] += (float) fabs(psi[FIX(i,j,k)]-var[LOCMIN][FIX(i,j,k)])/dens*eta/(rho*cp*(gx[FIX(i,j,k)]-gx[FIX(i-1,j,k)])*(gz[FIX(i,j,k)]-gz[FIX(i,j,k-1)])*(gy[FIX(i,j,k)]-gy[FIX(i,j-1,k)]));
      END_FOR
    }
    else {
      FOR_EACH_CELL
        if(flagp[FIX(i,j,k)]>0 ) continue;
          dens_s[FIX(i,j,k)]=(psi[FIX(i,j,k)]-var[LOCMAX][FIX(i,j,k)])/dens1*eta/(rho*cp*(gx[FIX(i,j,k)]-gx[FIX(i-1,j,k)])*(gz[FIX(i,j,k)]-gz[FIX(i,j,k-1)])*(gy[FIX(i,j,k)]-gy[FIX(i,j-1,k)]));
          psi[FIX(i,j,k)] += (float) fabs(psi[FIX(i,j,k)]-var[LOCMAX][FIX(i,j,k)])/dens1*eta/(rho*cp*(gx[FIX(i,j,k)]-gx[FIX(i-1,j,k)])*(gz[FIX(i,j,k)]-gz[FIX(i,j,k-1)])*(gy[FIX(i,j,k)]-gy[FIX(i,j-1,k)]));
      END_FOR
    }
  }
  else {
    mass0=0;

    FOR_EACH_CELL
      if(flagp[FIX(i,j,k)]<=0)  {
      var[LOCMIN][FIX(i,j,k)]=check_min(para, psi0, BINDEX[4][FIX(i,j,k)],  BINDEX[5][FIX(i,j,k)],  BINDEX[6][FIX(i,j,k)]); 
      var[LOCMAX][FIX(i,j,k)]=check_max(para, psi0, BINDEX[4][FIX(i,j,k)],  BINDEX[5][FIX(i,j,k)],  BINDEX[6][FIX(i,j,k)]); 
      }
    END_FOR

    qin= inflow(para,var,psi0,BINDEX);
    qout=outflow(para,var,psi0,BINDEX);
    mass0 +=(qin+qout)*dt;
    //printf("qin=%f\tqout=%f\n",qin,qout);
    //printf("mass0=%f\t%f\n",mass0,(qin+qout)*dt );

    FOR_EACH_CELL
      if(flagp[FIX(i,j,k)]>0) continue;
      dA= (gx[FIX(i,j,k)]-gx[FIX(i-1,j,k)])*(gz[FIX(i,j,k)]-gz[FIX(i,j,k-1)])*(gy[FIX(i,j,k)]-gy[FIX(i,j-1,k)]);
      area += dA;
      mass0 += rho*cp*psi0[FIX(i,j,k)]*dA;
      massstar += rho*cp*psi0[FIX(i,j,k)]*dA;
      mass  += rho*cp*psi[FIX(i,j,k)]*dA;
      dens +=  (float) fabs(psi[FIX(i,j,k)]-var[LOCMIN][FIX(i,j,k)]);
      dens1 +=  (float) fabs(psi[FIX(i,j,k)]-var[LOCMAX][FIX(i,j,k)]);   
    END_FOR

    massstar =mass0;
    eta=mass0-mass;
    //printf("eta=%f\n",eta);
    if(eta<0) {   // mass generation occured
      FOR_EACH_CELL
        if(flagp[FIX(i,j,k)]>0) continue;
        dens_s[FIX(i,j,k)]=(psi[FIX(i,j,k)]-var[LOCMIN][FIX(i,j,k)])/dens*eta/(rho*cp*(gx[FIX(i,j,k)]-gx[FIX(i-1,j,k)])*(gz[FIX(i,j,k)]-gz[FIX(i,j,k-1)])*(gy[FIX(i,j,k)]-gy[FIX(i,j-1,k)]));
        psi[FIX(i,j,k)] += (float) fabs(psi[FIX(i,j,k)]-var[LOCMIN][FIX(i,j,k)])/dens*eta/(rho*cp*(gx[FIX(i,j,k)]-gx[FIX(i-1,j,k)])*(gz[FIX(i,j,k)]-gz[FIX(i,j,k-1)])*(gy[FIX(i,j,k)]-gy[FIX(i,j-1,k)]));
      END_FOR
    }
    else {
      FOR_EACH_CELL
        if(flagp[FIX(i,j,k)]>0) continue;
          dens_s[FIX(i,j,k)]=(psi[FIX(i,j,k)]-var[LOCMAX][FIX(i,j,k)])/dens1*eta/(rho*cp*(gx[FIX(i,j,k)]-gx[FIX(i-1,j,k)])*(gz[FIX(i,j,k)]-gz[FIX(i,j,k-1)])*(gy[FIX(i,j,k)]-gy[FIX(i,j-1,k)]));
          psi[FIX(i,j,k)] += (float) fabs(psi[FIX(i,j,k)]-var[LOCMAX][FIX(i,j,k)])/dens1*eta/(rho*cp*(gx[FIX(i,j,k)]-gx[FIX(i-1,j,k)])*(gz[FIX(i,j,k)]-gz[FIX(i,j,k-1)])*(gy[FIX(i,j,k)]-gy[FIX(i,j-1,k)]));
      END_FOR
    }
  } 

  qin= inflow(para,var,psi,BINDEX);
  qout=outflow(para,var,psi,BINDEX);
  //printf("qin+qout=%f\n",(qin+qout) );
  mass=0;
  FOR_EACH_CELL
    if(flagp[FIX(i,j,k)]>0) continue;
    dA= (gx[FIX(i,j,k)]-gx[FIX(i-1,j,k)])*(gz[FIX(i,j,k)]-gz[FIX(i,j,k-1)])*(gy[FIX(i,j,k)]-gy[FIX(i,j-1,k)]);
    mass  += rho*cp*psi[FIX(i,j,k)]*dA;
  END_FOR
  //printf("mass-massstar=%f\n",mass-massstar);

} // End of mass_conservation()


REAL outflow(PARA_DATA *para, REAL **var, REAL *psi, int **BINDEX)
{
  int i, j, k;
  int it;
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
        //printf("Tout=%f\n",psi[FIX(SI,j+1,k+1)]);
    }
    if(SJ==EJ)
    {
        for(i=SI ;i<EI ;i++)
         for(k=SK ;k<EK ;k++)
         {
           mass_out += rho*cp*psi[FIX(i+1,SJ,k+1)]*flagp[FIX(i+1,SJ,k+1)]*v[FIX(i+1,SJ,k+1)]
                     *(gx[FIX(i+1,SJ,k+1)]-gx[FIX(i,SJ,k+1)])* (gz[FIX(i+1,SJ,k+1)]-gz[FIX(i+1,SJ,k)]);
         }
        //printf("Tout=%f\n",psi[FIX(i+1,SJ,k+1)]);
    }

        if(SK==EK)
    {
        for(i=SI ;i<EI ;i++)
          for(j=SJ ;j<EJ ;j++)
         {
           mass_out += rho*cp*psi[FIX(i+1,j+1,SK)]*flagp[FIX(i+1,j+1,SK)]*w[FIX(i+1,j+1,SK)]
                     *(gx[FIX(i+1,j+1,SK)]-gx[FIX(i,j+1,SK)])* (gy[FIX(i+1,j+1,SK)]-gy[FIX(i+1,j,SK)]);
           }
        //printf("Tout=%f\n",psi[FIX(i+1,j+1,SK)]);
    }
     
   }

   return mass_out;

}




REAL inflow(PARA_DATA *para, REAL **var, REAL *psi, int **BINDEX)
{
  int i, j, k;
  int it;
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


/***************inlet***************************************/

     for(zone_num=1;zone_num<=zone_inlet;zone_num++)
      {

     SI=BINDEX[7][zone_num];
     EI=BINDEX[8][zone_num];
     SJ=BINDEX[9][zone_num];
     EJ=BINDEX[10][zone_num];
     SK=BINDEX[11][zone_num];
     EK=BINDEX[12][zone_num];

        if(SI==EI)
    {
       for(j=SJ ;j<EJ ;j++)
         for(k=SK ;k<EK ;k++)
         {
              
           mass_in +=rho*cp*para->bc->t_bc[zone_num]*flagp[FIX(SI,j+1,k+1)]*u[FIX(SI,j+1,k+1)]
                     *(gy[FIX(SI,j+1,k+1)]-gy[FIX(SI,j,k+1)])* (gz[FIX(SI,j+1,k+1)]-gz[FIX(SI,j+1,k)]);
         }
    }
        if(SJ==EJ)
    {
        for(i=SI ;i<EI ;i++)
         for(k=SK ;k<EK ;k++)
         {
           mass_in +=rho*cp*para->bc->t_bc[zone_num]* flagp[FIX(i+1,SJ,k+1)]*v[FIX(i+1,SJ,k+1)]
                     *(gx[FIX(i+1,SJ,k+1)]-gx[FIX(i,SJ,k+1)])* (gz[FIX(i+1,SJ,k+1)]-gz[FIX(i+1,SJ,k)]);
         }
    }

        if(SK==EK)
    {
        for(i=SI ;i<EI ;i++)
          for(j=SJ ;j<EJ ;j++)
         {
           mass_in += rho*cp*para->bc->t_bc[zone_num]*flagp[FIX(i+1,j+1,SK)]*w[FIX(i+1,j+1,SK)]
                     *(gx[FIX(i+1,j+1,SK)]-gx[FIX(i,j+1,SK)])* (gy[FIX(i+1,j+1,SK)]-gy[FIX(i+1,j,SK)]);
           }
    }
    
   }
  return mass_in;
}


REAL check_min(PARA_DATA *para, REAL *psi, int ci,int cj,int ck)
{
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

}// End of check_residual( )


REAL check_max( PARA_DATA *para, REAL *psi, int ci,int cj,int ck)
{
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

}// End of check_residual( )

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

}// End of check_residual( )


REAL check_max_pix( PARA_DATA *para, REAL *psi)
{
  int i, j, k;
  REAL tmp=psi[PIX(0,0,0)];

  for(i=0;i<=1;i++)
    for(j=0;j<=1;j++)
          for(k=0;k<=1;k++)
        {
          if(tmp<psi[PIX(i,j,k)]) tmp=psi[PIX(i,j,k)];
        
        }
return tmp;

}// End of check_residual( )


REAL check_avg(PARA_DATA *para, REAL **var, REAL *psi)
{
  int imax = para->geom->imax, jmax = para->geom->jmax; 
  int kmax = para->geom->kmax;
  int i, j, k;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  REAL tmp=0,vol=0;
  REAL *flagp=var[FLAGP];
  REAL *gx=var[GX],*gy=var[GY],*gz=var[GZ];
  REAL dx,dy,dz;
  REAL rho=para->prob->rho;
  REAL spec=para->prob->spec;
  

  FOR_EACH_CELL
    if(flagp[FIX(i,j,k)]==1) continue;

    dx=gx[FIX(i,j,k)]-gx[FIX(i-1,j,k)];
    dy=gy[FIX(i,j,k)]-gy[FIX(i,j-1,k)];
    dz=gz[FIX(i,j,k)]-gz[FIX(i,j,k-1)];
     tmp += dx*dy*dz*(float) fabs(psi[FIX(i,j,k)]);

   vol += dx*dy*dz;
  END_FOR
    
  //return tmp / vol;
  return tmp*rho*spec;

}// End of check_residual( )




REAL qwall(PARA_DATA *para, REAL **var,int **BINDEX)
{
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

  for(it=1;it<=indexmax;it++)
      {
        i=BINDEX[0][it];
        j=BINDEX[1][it];
        k=BINDEX[2][it];
    
      /*******X direction EAST FACE********/

        if(i<imax && flagp[FIX(i+1,j,k)]<0)
           {
           switch((int)flagu[FIX(i,j,k)])
           {
              case 0: //inlet_surface
            case 1: //outlet_surface
              qwall_temp=0;
              break;
            case 2: //wall_surface
                         zone_num=BINDEX[13][FIX(i,j,k)];
              if(para->bc->fltmp[zone_num]==0)
               {
               qflux=para->bc->t_bc[zone_num];
               qwall_temp = qflux*(gy[FIX(i,j,k)]-gy[FIX(i,j-1,k)])*(gz[FIX(i,j,k)]-gz[FIX(i,j,k-1)]);
               qwall += qwall_temp;
       
               }
              else
               {
                qwall_temp=rho*cp*(para->bc->t_bc[zone_num]-psi[FIX(i+1,j,k)])*coeff_h*(gy[FIX(i,j,k)]-gy[FIX(i,j-1,k)])*(gz[FIX(i,j,k)]-gz[FIX(i,j,k-1)]);
                        //  if(qwall_temp<0) printf("zonenum=%d\ttblock=%f\n",BINDEX[16][FIX(i,j,k)],para->bc->t_bc[zone_num]-psi[FIX(i+1,j,k)]);
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
               qwall_temp = rho*cp*(para->bc->t_bc[zone_num]-psi[FIX(i+1,j,k)])*coeff_h*(gy[FIX(i,j,k)]-gy[FIX(i,j-1,k)])*(gz[FIX(i,j,k)]-gz[FIX(i,j,k-1)]);    
              // if(qwall_temp<0) printf("zonenum=%d\ttblock=%f\t%f\t%f\t%d\n",BINDEX[16][FIX(i,j,k)],para->bc->t_bc[zone_num],gy[FIX(i,j,k)]-gy[FIX(i,j-1,k)],gz[FIX(i,j,k)]-gz[FIX(i,j,k-1)],k);
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
               qwall_temp = rho*cp*(para->bc->t_bc[zone_num]-psi[FIX(i-1,j,k)])*coeff_h*(gy[FIX(i,j,k)]-gy[FIX(i,j-1,k)])*(gz[FIX(i,j,k)]-gz[FIX(i,j,k-1)]);
                      //    if(qwall_temp<0) printf("zonenum=%d\ttblock=%f\n",BINDEX[16][FIX(i,j,k)],para->bc->t_bc[zone_num]-psi[FIX(i-1,j,k)]); 
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
               qwall_temp = rho*cp*(para->bc->t_bc[zone_num]-psi[FIX(i-1,j,k)])*coeff_h*(gy[FIX(i,j,k)]-gy[FIX(i,j-1,k)])*(gz[FIX(i,j,k)]-gz[FIX(i,j,k-1)]);
              // if(qwall_temp<0) printf("zonenum=%d\ttblock=%f\t%f\t%f\t%d\n",BINDEX[16][FIX(i,j,k)],para->bc->t_bc[zone_num],gy[FIX(i,j,k)]-gy[FIX(i,j-1,k)],gz[FIX(i,j,k)]-gz[FIX(i,j,k-1)],k);
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
               qwall_temp = rho*cp*(para->bc->t_bc[zone_num]-psi[FIX(i,j+1,k)])*coeff_h*(gx[FIX(i,j,k)]-gx[FIX(i-1,j,k)])*(gz[FIX(i,j,k)]-gz[FIX(i,j,k-1)]);
                       //   if(qwall_temp<0) printf("zonenum=%d\ttblock=%f\n",BINDEX[16][FIX(i,j,k)],para->bc->t_bc[zone_num]-psi[FIX(i,j+1,k)]); 
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
               qwall_temp = rho*cp*(para->bc->t_bc[zone_num]-psi[FIX(i,j+1,k)])*coeff_h*(gx[FIX(i,j,k)]-gx[FIX(i-1,j,k)])*(gz[FIX(i,j,k)]-gz[FIX(i,j,k-1)]);
              // if(qwall_temp<0) printf("zonenum=%d\ttblock=%f%d\t%d\t%d\n",BINDEX[16][FIX(i,j,k)],para->bc->t_bc[zone_num],i,j,k);
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
                 qwall_temp = rho*cp*(para->bc->t_bc[zone_num]-psi[FIX(i,j-1,k)])*coeff_h*(gx[FIX(i,j,k)]-gx[FIX(i-1,j,k)])*(gz[FIX(i,j,k)]-gz[FIX(i,j,k-1)]);
                     //     if(qwall_temp<0) printf("zonenum=%d\ttblock=%f\n",BINDEX[16][FIX(i,j,k)],para->bc->t_bc[zone_num]-psi[FIX(i,j-1,k)]); 
               qwall += qwall_temp;           
               }
              break;
            case 3:
                           zone_num=BINDEX[16][FIX(i,j,k)];
                // printf("zone_num=%d temperatue=%f\n",zone_num, para->bc->t_bc[zone_num]);
               if(para->bc->fltmp[zone_num]==0)
               {
               qflux=para->bc->t_bc[zone_num];
               qwall_temp = qflux*(gx[FIX(i,j,k)]-gx[FIX(i-1,j,k)])*(gz[FIX(i,j,k)]-gz[FIX(i,j,k-1)]); 
               qwall += qwall_temp;
               }
              else
               {
                qwall_temp = rho*cp*(para->bc->t_bc[zone_num]-psi[FIX(i,j-1,k)])*coeff_h*(gx[FIX(i,j,k)]-gx[FIX(i-1,j,k)])*(gz[FIX(i,j,k)]-gz[FIX(i,j,k-1)]);
              // if(qwall_temp<0) printf("zonenum=%d\ttblock=%f%d\t%d\t%d\n",BINDEX[16][FIX(i,j,k)],para->bc->t_bc[zone_num],i,j,k);
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
                qwall_temp = rho*cp*(para->bc->t_bc[zone_num]-psi[FIX(i,j,k+1)])*coeff_h*(gy[FIX(i,j,k)]-gy[FIX(i,j-1,k)])*(gx[FIX(i,j,k)]-gx[FIX(i-1,j,k)]);
                       //   if(qwall_temp<0) printf("zonenum=%d\ttblock=%f\n",BINDEX[16][FIX(i,j,k)],para->bc->t_bc[zone_num]-psi[FIX(i,j,k+1)]); 
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
                qwall_temp = rho*cp*(para->bc->t_bc[zone_num]-psi[FIX(i,j,k+1)])*coeff_h*(gy[FIX(i,j,k)]-gy[FIX(i,j-1,k)])*(gx[FIX(i,j,k)]-gx[FIX(i-1,j,k)]);
            //  if(qwall_temp<0) printf("zonenum=%d\ttblock=%f\t%f\t%f\t%f\n",BINDEX[16][FIX(i,j,k)],para->bc->t_bc[zone_num]-psi[FIX(i,j,k+1)],gy[FIX(i,j,k)]-gy[FIX(i,j-1,k)],gx[FIX(i,j,k)]-gx[FIX(i-1,j,k)],j);
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
               qwall_temp = rho*cp*(para->bc->t_bc[zone_num]-psi[FIX(i,j,k-1)])*coeff_h*(gy[FIX(i,j,k)]-gy[FIX(i,j-1,k)])*(gx[FIX(i,j,k)]-gx[FIX(i-1,j,k)]);
                        //  if(qwall_temp<0) printf("zonenum=%d\ttblock=%f\n",BINDEX[16][FIX(i,j,k)],para->bc->t_bc[zone_num]-psi[FIX(i,j,k-1)]);  
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
                qwall_temp = rho*cp*(para->bc->t_bc[zone_num]-psi[FIX(i,j,k-1)])*coeff_h*(gy[FIX(i,j,k)]-gy[FIX(i,j-1,k)])*(gx[FIX(i,j,k)]-gx[FIX(i-1,j,k)]);
              // if(qwall_temp<0) printf("zonenum=%d\ttblock=%f%d\t%d\t%d\n",BINDEX[16][FIX(i,j,k)],para->bc->t_bc[zone_num],i,j,k);
               qwall += qwall_temp;
          
               }
               break;
           }
           qb +=qwall_temp;
      }

    }//end_for

         //printf("qwall=%f\n",qwall);

      /**********************************Heat sources****************************************************/
      //      for(it=indexmax+1;it<=indexmax_us;it++)
      //{
      //  i=BINDEX[0][it];
      //          j=BINDEX[1][it];
      //          k=BINDEX[2][it];
      //  zone_num=BINDEX[16][FIX(i,j,k)]; 
      //  qwall += para->bc->qs;   //para->bc->t_bc[zone_num];

      //}
       qwall += para->bc->qs;   //para->bc->t_bc[zone_num];
  //printf("qwall=%f\t%f\t%f\t%f\t%f\t%f\n",qw,qe,qs,qn,qb,qf);
  
  return qwall;
      

}

REAL check_energy(PARA_DATA *para, REAL **var,int **BINDEX)
{
  REAL qin,qout;
  REAL qw=0.0000001f;
    REAL error;

  qin= inflow(para, var, var[TEMP], BINDEX);
  qout=outflow(para, var, var[TEMP], BINDEX);
  qw += qwall(para, var, BINDEX);

  error=(qout+qin+qw)/qw*100;

 // printf("qin=%f\tqout=%f\tqw=%f\terror=%f\n",qin,qout,qw,error);

return error;  
  

}

void swapuvw(PARA_DATA *para, REAL **var)
{
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

}


void swapuvw_temp(PARA_DATA *para, REAL **var)
{
int i, j, k;
  int imax = para->geom->imax, jmax = para->geom->jmax;
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
   REAL *u1 = var[VXC], *v1 = var[VYC], *w1 = var[VZC];  
  REAL *u = var[VX], *v = var[VY], *w = var[VZ];  
  

  FOR_ALL_CELL
    u1[FIX(i,j,k)]=u[FIX(i,j,k)];
    v1[FIX(i,j,k)]=v[FIX(i,j,k)];
    w1[FIX(i,j,k)]=w[FIX(i,j,k)];
  END_FOR

}

REAL num_visco(REAL dx, REAL dt,REAL x_1)
{

 return x_1*(1-x_1)*dx*dx/(2.0f*dt);

} 

void cfl(PARA_DATA *para, REAL **var)
{
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

}


REAL velmag(REAL u, REAL v)
{
 return sqrt(u*u+v*v);
}

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
           u[FIX(i,j,k)]= interpolation_linear(y_1, z_1, u0[FIX(i,j,k)], u0[FIX(i,j,k+1)], u0[FIX(i,j+1,k)], u0[FIX(i,j+1,k+1)]);
        
        }
        }

}


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
           v[FIX(i,j,k)]= interpolation_linear(y_1, z_1, v0[FIX(i,j,k)], v0[FIX(i,j,k+1)], v0[FIX(i+1,j,k)], v0[FIX(i+1,j,k+1)]);
        
        }
        }
}

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
           w[FIX(i,j,k)]= interpolation_linear(y_1, z_1, w0[FIX(i,j,k)], w0[FIX(i+1,j,k)], w0[FIX(i,j+1,k)], w0[FIX(i+1,j+1,k)]);
        
        }
        }

}


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
                                             p0[FIX(i,j,k)],  p0[FIX(i,j+1,k)],  p0[FIX(i+1,j,k)],  p0[FIX(i+1,j+1,k)], 
                                             p0[FIX(i,j,k+1)],p0[FIX(i,j+1,k+1)],p0[FIX(i+1,j,k+1)],p0[FIX(i+1,j+1,k+1)]);
        
        }
        

    i=0;
    for(j=1;j<jmax;j++)
       for(k=1;k<kmax;k++)
       {
           y_1= (gy[FIX(i,j,k)]-y[FIX(i,j,k)])/(y[FIX(i,j+1,k)]-y[FIX(i,j,k)]);
                   z_1= (gz[FIX(i,j,k)]-z[FIX(i,j,k)])/(z[FIX(i,j,k+1)]-z[FIX(i,j,k)]);
           p[FIX(i,j,k)]= interpolation_linear(y_1, z_1, p0[FIX(i,j,k)], p0[FIX(i,j,k+1)], p0[FIX(i,j+1,k)], p0[FIX(i,j+1,k+1)]);
              
         }

    i=imax;
    for(j=1;j<jmax;j++)
       for(k=1;k<kmax;k++)
       {
           y_1= (gy[FIX(i,j,k)]-y[FIX(i,j,k)])/(y[FIX(i,j+1,k)]-y[FIX(i,j,k)]);
                   z_1= (gz[FIX(i,j,k)]-z[FIX(i,j,k)])/(z[FIX(i,j,k+1)]-z[FIX(i,j,k)]);
           p[FIX(i,j,k)]= interpolation_linear(y_1, z_1, p0[FIX(i,j,k)], p0[FIX(i,j,k+1)], p0[FIX(i,j+1,k)], p0[FIX(i,j+1,k+1)]);
              
         }

  j=0;
      for(i=1;i<imax;i++)
       for(k=1;k<kmax;k++)
       {
             y_1= (gx[FIX(i,j,k)]-x[FIX(i,j,k)])/(x[FIX(i+1,j,k)]-x[FIX(i,j,k)]);
                   z_1= (gz[FIX(i,j,k)]-z[FIX(i,j,k)])/(z[FIX(i,j,k+1)]-z[FIX(i,j,k)]);
           p[FIX(i,j,k)]= interpolation_linear(y_1, z_1, p0[FIX(i,j,k)], p0[FIX(i,j,k+1)], p0[FIX(i+1,j,k)], p0[FIX(i+1,j,k+1)]);
       }
  j=jmax;
      for(i=1;i<imax;i++)
       for(k=1;k<kmax;k++)
       {
             y_1= (gx[FIX(i,j,k)]-x[FIX(i,j,k)])/(x[FIX(i+1,j,k)]-x[FIX(i,j,k)]);
                   z_1= (gz[FIX(i,j,k)]-z[FIX(i,j,k)])/(z[FIX(i,j,k+1)]-z[FIX(i,j,k)]);
           p[FIX(i,j,k)]= interpolation_linear(y_1, z_1, p0[FIX(i,j,k)], p0[FIX(i,j,k+1)], p0[FIX(i+1,j,k)], p0[FIX(i+1,j,k+1)]);
       }
     k=0;
       for(i=1;i<imax;i++)
          for(j=1;j<jmax;j++)
        {
           y_1= (gy[FIX(i,j,k)]-y[FIX(i,j,k)])/(y[FIX(i,j+1,k)]-y[FIX(i,j,k)]);
                   z_1= (gx[FIX(i,j,k)]-x[FIX(i,j,k)])/(x[FIX(i+1,j,k)]-x[FIX(i,j,k)]);
           p[FIX(i,j,k)]= interpolation_linear(y_1, z_1, p0[FIX(i,j,k)], p0[FIX(i+1,j,k)], p0[FIX(i,j+1,k)], p0[FIX(i+1,j+1,k)]);
        
        }
     k=kmax;
       for(i=1;i<imax;i++)
          for(j=1;j<jmax;j++)
        {
           y_1= (gy[FIX(i,j,k)]-y[FIX(i,j,k)])/(y[FIX(i,j+1,k)]-y[FIX(i,j,k)]);
                   z_1= (gx[FIX(i,j,k)]-x[FIX(i,j,k)])/(x[FIX(i+1,j,k)]-x[FIX(i,j,k)]);
           p[FIX(i,j,k)]= interpolation_linear(y_1, z_1, p0[FIX(i,j,k)], p0[FIX(i+1,j,k)], p0[FIX(i,j+1,k)], p0[FIX(i+1,j+1,k)]);
        
        }

}

void corners(PARA_DATA *para, REAL **var, REAL *psi)
{
  int i,j,k;
  int imax=para->geom->imax, jmax=para->geom->jmax;
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);



   for(k=1;k<=kmax;k++)
    {
    psi[FIX(0,0,k)]   = (psi[FIX(1,0,k)]+psi[FIX(0,1,k)])/2.0f;
    psi[FIX(0,jmax+1,k)]= (psi[FIX(1,jmax+1,k)]+psi[FIX(0,jmax,k)])/2.0f;
    psi[FIX(imax+1,0,k)]= (psi[FIX(imax,0,k)]+psi[FIX(imax+1,1,k)])/2.0f;
    psi[FIX(imax+1,jmax+1,k)]=(psi[FIX(imax,jmax+1,k)]+psi[FIX(imax+1,jmax,k)])/2.0f;
    }

  for(j=1;j<=jmax;j++)
    {
    psi[FIX(0,j,0)]   = (psi[FIX(1,j,0)]+psi[FIX(0,j,1)])/2.0f;
    psi[FIX(0,j,kmax+1)]= (psi[FIX(1,j,kmax+1)]+psi[FIX(0,j,kmax)])/2.0f;
    psi[FIX(imax+1,j,0)]= (psi[FIX(imax,j,0)]+psi[FIX(imax+1,j,1)])/2.0f;
    psi[FIX(imax+1,j,kmax+1)]=(psi[FIX(imax,j,kmax+1)]+psi[FIX(imax+1,j,kmax)])/2.0f;
    }

    for(i=1;i<=imax;i++)
    {
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


} //corners()

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


} //corners()



void uvw_center(PARA_DATA *para, REAL **var)
{
  int i, j, k;
  int imax=para->geom->imax, jmax=para->geom->jmax;
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  REAL *u0 = var[VX], *v0 = var[VY], *w0 = var[VZ];
  REAL *flagp=var[FLAGP];


      for(j=0; j<=jmax+1; j++)
    {
    for(k=0; k<=kmax+1; k++)
          {
            
                u0[FIX(imax+1,j,k)] = u0[FIX(imax,j,k)];

                     for(i=imax; i>=1; i--)
                     {
                         if(flagp[FIX(i,j,k)]==1) u0[FIX(i,j,k)]=0;
             else u0[FIX(i,j,k)] = 0.5f * (u0[FIX(i,j,k)]+u0[FIX(i-1,j,k)]);
                      }
            }
    }


  
    for(i=0; i<=imax+1; i++)
    {
    for(k=0; k<=kmax+1; k++)
          {
            
                v0[FIX(i,jmax+1,k)] = v0[FIX(i,jmax,k)];

                     for(j=jmax; j>=1; j--)
                     {
                         if(flagp[FIX(i,j,k)]==1) v0[FIX(i,j,k)]=0;
             else v0[FIX(i,j,k)] = 0.5f * (v0[FIX(i,j,k)]+v0[FIX(i,j-1,k)]);
                        }
            }
    }
   

   for(i=0; i<=imax+1; i++)
    {
    for(j=0; j<=jmax+1; j++)
          {
            
                w0[FIX(i,j,kmax+1)] = w0[FIX(i,j,kmax)];

                     for(k=kmax; k>=1; k--)
                     {
              if(flagp[FIX(i,j,k)]==1) w0[FIX(i,j,k)]=0;
                          else w0[FIX(i,j,k)] = 0.5f * (w0[FIX(i,j,k)]+w0[FIX(i,j,k-1)]);
                      }
            }
    }
}


REAL qheatsource(PARA_DATA *para, REAL **var,int **BINDEX)
{
  int i, j, k;
  int it;
  int indexmax,indexmax_us,zone_num;
  int imax = para->geom->imax, jmax = para->geom->jmax;
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  REAL *T=var[TEMP];
  REAL *vt=var[VT];
  REAL *gx = var[GX], *gy = var[GY], *gz = var[GZ];
  REAL *x = var[X], *y = var[Y], *z = var[Z];
  REAL qflux=0;
  REAL dxe, dxw, dyn, dys, dzf, dzb, Dx, Dy, Dz;
  REAL ae, aw, an,as,af,ab;
  REAL kapa = para->prob->nu; 
  REAL rho=para->prob->rho;
  REAL cp=para->prob->spec;
  REAL Prt= para->prob->Prt;
 
 indexmax=para->geom->index[1];
 indexmax_us=para->geom->index[2];
 i=BINDEX[0][indexmax_us];
 k=BINDEX[2][indexmax_us];

 if(para->prob->plume_mod==1) {

 for(j=BINDEX[1][indexmax_us];j<(BINDEX[1][indexmax_us]+para->geom->plmax);j++) {
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

          aw= rho*cp*(kapa+0.5f*(vt[FIX(i,j,k)]+vt[FIX(i-1,j,k)]))/Prt*Dy*Dz/dxw;
          ae= rho*cp*(kapa+0.5f*(vt[FIX(i,j,k)]+vt[FIX(i+1,j,k)]))/Prt*Dy*Dz/dxe;
          an= rho*cp*(kapa+0.5f*(vt[FIX(i,j,k)]+vt[FIX(i,j+1,k)]))/Prt*Dx*Dz/dyn;
          as= rho*cp*(kapa+0.5f*(vt[FIX(i,j,k)]+vt[FIX(i,j-1,k)]))/Prt*Dx*Dz/dys;
          af= rho*cp*(kapa+0.5f*(vt[FIX(i,j,k)]+vt[FIX(i,j,k+1)]))/Prt*Dx*Dy/dzf;
          ab= rho*cp*(kapa+0.5f*(vt[FIX(i,j,k)]+vt[FIX(i,j,k-1)]))/Prt*Dx*Dy/dzb;
      /*   printf("aw=%f\t%f\t%f\t%f\t%f\t%f\n", aw,ae,an,as,af,ab);
         printf("Tw=%f\t%f\t%f\t%f\t%f\t%f\n", T[FIX(i,j,k)]-T[FIX(i-1,j,k)],
                                               T[FIX(i,j,k)]-T[FIX(i+1,j,k)],
                                               T[FIX(i,j,k)]-T[FIX(i,j+1,k)],                                                
                                               T[FIX(i,j,k)]-T[FIX(i,j-1,k)],                                              
                                               T[FIX(i,j,k)]-T[FIX(i,j,k+1)], 
                                               T[FIX(i,j,k)]-T[FIX(i,j,k-1)]);
                                               */
        }
        else {
          aw= rho*cp*kapa/Prt*Dy*Dz/dxw;
          ae= rho*cp*kapa/Prt*Dy*Dz/dxe;
          an= rho*cp*kapa/Prt*Dx*Dz/dyn;
          as= rho*cp*kapa/Prt*Dx*Dz/dys;
          af= rho*cp*kapa/Prt*Dx*Dy/dzf;
          ab= rho*cp*kapa/Prt*Dx*Dy/dzb;
        }

        
          aw= 0;
          ae= 0;
          af= 0;
          ab= 0;
          an=0;
          as=0;


        //qflux=para->bc->qflow_diff;

        qflux +=( ae*(T[FIX(i,j,k)]-T[FIX(i+1,j,k)])
              +aw*(T[FIX(i,j,k)]-T[FIX(i-1,j,k)])
              +af*(T[FIX(i,j,k)]-T[FIX(i,j,k+1)])
              +ab*(T[FIX(i,j,k)]-T[FIX(i,j,k-1)]))
              +an*(T[FIX(i,j,k)]-T[FIX(i,j+1,k)])
              +as*(T[FIX(i,j,k)]-T[FIX(i,j-1,k)]);


  }
  //printf("qdiff=%f\t%f\n",para->bc->qflow_diff,qflux);
 }
 else {
   qflux=0;
 }

  return qflux;
}