#include <stdio.h>
#include <stdlib.h>
#include "data_structure.h"
#include "boundary.h"
#include "cprojection.h"
#include "projection.h"
#include "solver_gs.h"
#include "utility.h"
#include "interpolation.h"
#include <math.h>

void cprojection(PARA_DATA *para, REAL **var, int clvl, int **BINDEX)
{
  int i,j,k;
  int NMG=para->geom->NMG;
  int lvl;
  int mglvl;
  int NI[10],NJ[10],NK[10],IMAX[10],IJMAX[10],IJK[10];
  REAL *x = var[X], *y = var[Y], *z = var[Z];
  REAL *u=var[VX],*v=var[VY],*w=var[VZ];
  REAL *gx = var[GX], *gy = var[GY], *gz = var[GZ]; 
  REAL *b = var[B], *ap = var[AP], *ab = var[AB], *af = var[AF];
  REAL *ae = var[AE], *aw =var[AW], *an = var[AN], *as = var[AS];
  REAL residual = 1.0;  
  REAL dt=para->mytime->dt;
  REAL dxe,dxw, dyn,dys,dzf,dzb,Dx,Dy,Dz;

   mglvl=NMG-clvl;

   NI[mglvl]=para->geom->NI[mglvl];
   NJ[mglvl]=para->geom->NJ[mglvl]; 
   NK[mglvl]=para->geom->NK[mglvl]; 
   IMAX[mglvl]=para->geom->IMAX[mglvl];
   IJMAX[mglvl]=para->geom->IJMAX[mglvl];
   IJK[mglvl]=para->geom->IJK[mglvl];



  for(lvl=NMG-1;lvl>=mglvl;lvl--) cmap (para, var,lvl);

  for(i=1;i<=NI[mglvl];i++)
      for(j=1;j<=NJ[mglvl];j++)
        for(k=1;k<=NK[mglvl];k++)
      {
                    dxe= x[MIX(i+1,j,k,mglvl)]-x[MIX(i,j,k,mglvl)];
          dxw= x[MIX(i,j,k,mglvl)]-x[MIX(i-1,j,k,mglvl)];
          dyn= y[MIX(i,j+1,k,mglvl)]-y[MIX(i,j,k,mglvl)];
          dys= y[MIX(i,j ,k,mglvl)]-y[MIX(i,j-1,k,mglvl)];
          dzf= z[MIX(i,j ,k+1,mglvl)]-z[MIX(i,j,k,mglvl)];
          dzb= z[MIX(i,j ,k,mglvl)]-z[MIX(i,j,k-1,mglvl)];
          Dx  = gx[MIX(i,j,k,mglvl)]-gx[MIX(i-1,j,k,mglvl)];
          Dy  = gy[MIX(i,j ,k,mglvl)]-gy[MIX(i,j-1,k,mglvl)];
          Dz  = gz[MIX(i,j ,k,mglvl)]-gz[MIX(i,j,k-1,mglvl)];

 
           ae[MIX(i,j,k,mglvl)] = Dy*Dz/dxe;
           aw[MIX(i,j,k,mglvl)] = Dy*Dz/dxw;      
          an[MIX(i,j,k,mglvl)] = Dx*Dz/dyn;
          as[MIX(i,j,k,mglvl)] = Dx*Dz/dys;
          af[MIX(i,j,k,mglvl)] = Dx*Dy/dzf;
          ab[MIX(i,j,k,mglvl)] = Dx*Dy/dzb;
          b[MIX(i,j,k,mglvl)] = Dx*Dy*Dz/dt*((u[MIX(i-1,j,k,mglvl)]-u[MIX(i,j,k,mglvl)])/Dx
                    + (v[MIX(i,j-1,k,mglvl)]-v[MIX(i,j,k,mglvl)])/Dy
                      + (w[MIX(i,j,k-1,mglvl)]-w[MIX(i,j,k,mglvl)])/Dz);

      }

 set_bnd_pressure_mg(para, var, var[IP],mglvl,BINDEX);

  for(i=1;i<=NI[mglvl];i++)
      for(j=1;j<=NJ[mglvl];j++)
        for(k=1;k<=NK[mglvl];k++)
      {
                       ap[MIX(i,j,k,mglvl)] = ae[MIX(i,j,k,mglvl)] + aw[MIX(i,j,k,mglvl)] + as[MIX(i,j,k,mglvl)] + an[MIX(i,j,k,mglvl)]
                               + af[MIX(i,j,k,mglvl)] + ab[MIX(i,j,k,mglvl)];
      }

    GS_PMG(para, var, var[IP],mglvl);
  // SOR_PMG(para, var, var[IP],mglvl);
 
  //set_bnd_pressurevalue_mg(para, var, var[IP],mglvl,BINDEX);

 //for(lvl=mglvl+1;lvl<=NMG;lvl++) fmap (para, var,lvl);

// vel_correction( para,  var);

   vel_map (para,var,NMG);

} // End of projecttion( )


void fmap (PARA_DATA *para, REAL **var, int mglvl)
{

    int i,j,k;
  int p,q,r;
  int itx,ity,itz;
  int imax,jmax,kmax;
  int IND;
    int NMG=para->geom->NMG;
    int lvl;
    int NI[10],NJ[10],NK[10],IMAX[10],IJMAX[10],IJK[10];
  REAL *fx=var[FX],*fy=var[FY],*fz=var[FZ];
  REAL *flagp=var[FLAGP];



  for(lvl=1;lvl<=NMG;lvl++)
  {
   NI[lvl]=para->geom->NI[lvl];
   NJ[lvl]=para->geom->NJ[lvl]; 
   NK[lvl]=para->geom->NK[lvl]; 
   IMAX[lvl]=para->geom->IMAX[lvl];
   IJMAX[lvl]=para->geom->IJMAX[lvl];
   IJK[lvl]=para->geom->IJK[lvl];
  }

  /* corner     */

  imax=NI[mglvl-1];
  jmax=NJ[mglvl-1];
  kmax=NK[mglvl-1];

   for(i=1;i<=imax;i++)
  {
     var[IP][MIX(i,0,0, mglvl-1)]=(var[IP][MIX(i,0,1,mglvl-1)]+var[IP][MIX(i,1,0,mglvl-1)])/2.0f;
     var[IP][MIX(i,0,kmax, mglvl-1)]=(var[IP][MIX(i,0,kmax-1,mglvl-1)]+var[IP][MIX(i,1,kmax,mglvl-1)])/2.0f;  
     var[IP][MIX(i,jmax,0, mglvl-1)]=(var[IP][MIX(i,jmax,1,mglvl-1)]+var[IP][MIX(i,jmax+1,0,mglvl-1)])/2.0f;  
     var[IP][MIX(i,jmax,kmax, mglvl-1)]=(var[IP][MIX(i,jmax,kmax+1,mglvl-1)]+var[IP][MIX(i,jmax+1,kmax,mglvl-1)])/2.0f;  
  }
    for(j=1;j<=jmax;j++)
  {
     var[IP][MIX(0,j,0, mglvl-1)]=(var[IP][MIX(0,j,1,mglvl-1)]+var[IP][MIX(1,j,0,mglvl-1)])/2.0f;
     var[IP][MIX(imax,j,0, mglvl-1)]=(var[IP][MIX(imax,j,1,mglvl-1)]+var[IP][MIX(imax+1,j,0,mglvl-1)])/2.0f;
     var[IP][MIX(0,j,kmax, mglvl-1)]=(var[IP][MIX(0,j,kmax+1,mglvl-1)]+var[IP][MIX(1,j,kmax,mglvl-1)])/2.0f;
     var[IP][MIX(imax,j,kmax, mglvl-1)]=(var[IP][MIX(imax,j,kmax+1,mglvl-1)]+var[IP][MIX(imax+1,j,kmax,mglvl-1)])/2.0f;
  }
      for(k=1;k<=kmax;k++)
  {
     var[IP][MIX(0,0,k, mglvl-1)]=(var[IP][MIX(0,1,k,mglvl-1)]+var[IP][MIX(1,0,k,mglvl-1)])/2.0f;
     var[IP][MIX(imax,0,k, mglvl-1)]=(var[IP][MIX(imax,1,k,mglvl-1)]+var[IP][MIX(imax+1,0,k,mglvl-1)])/2.0f;
     var[IP][MIX(0,jmax,k, mglvl-1)]=(var[IP][MIX(0,jmax+1,k,mglvl-1)]+var[IP][MIX(1,jmax,k,mglvl-1)])/2.0f;
     var[IP][MIX(imax,jmax,k, mglvl-1)]=(var[IP][MIX(imax,jmax+1,k,mglvl-1)]+var[IP][MIX(imax+1,jmax,k,mglvl-1)])/2.0f;
  }


  var[IP][MIX(0,0,0, mglvl-1)] = (var[IP][MIX(0,1,0, mglvl-1)]+var[IP][MIX(1,0,0, mglvl-1)]+var[IP][MIX(0,0,1, mglvl-1)]) / 3.0f;

  var[IP][MIX(0,jmax+1,0, mglvl-1)] = ( var[IP][MIX(1,jmax+1,0, mglvl-1)]+var[IP][MIX(0,jmax,0, mglvl-1)]  +var[IP][MIX(0,jmax+1,1, mglvl-1)]) / 3.0f;

  var[IP][MIX(imax+1,0,0, mglvl-1)] = ( var[IP][MIX(imax,0,0, mglvl-1)]+var[IP][MIX(imax+1,1,0, mglvl-1)]  +var[IP][MIX(imax+1,0,1, mglvl-1)]) / 3.0f;

  var[IP][MIX(imax+1,jmax+1,0, mglvl-1)] = ( var[IP][MIX(imax,jmax+1,0, mglvl-1)]+var[IP][MIX(imax+1,jmax,0, mglvl-1)] +var[IP][MIX(imax+1,jmax+1,1, mglvl-1)]) / 3.0f;

  var[IP][MIX(0,0,kmax+1, mglvl-1)] = ( var[IP][MIX(0,1,kmax+1, mglvl-1)]+var[IP][MIX(1,0,kmax+1, mglvl-1)] +var[IP][MIX(0,0,kmax, mglvl-1)]) / 3.0f;  

  var[IP][MIX(0,jmax+1,kmax+1, mglvl-1)] = ( var[IP][MIX(1,jmax+1,kmax+1, mglvl-1)]+var[IP][MIX(0,jmax,kmax+1, mglvl-1)]+var[IP][MIX(0,jmax+1,kmax, mglvl-1)]) / 3.0f;

  var[IP][MIX(imax+1,0,kmax+1, mglvl-1)] = ( var[IP][MIX(imax,0,kmax+1, mglvl-1)]+var[IP][MIX(imax+1,1,kmax+1, mglvl-1)] +var[IP][MIX(imax+1,0,kmax, mglvl-1)]) / 3.0f;

  var[IP][MIX(imax+1,jmax+1,kmax+1, mglvl-1)] = ( var[IP][MIX(imax,jmax+1,kmax+1, mglvl-1)]+var[IP][MIX(imax+1,jmax,kmax+1, mglvl-1)]+var[IP][MIX(imax+1,jmax+1,kmax, mglvl-1)]) / 3.0f;
  /* corner     */


  for(i=1;i<=NI[mglvl-1];i++)
      for(j=1;j<=NJ[mglvl-1];j++)
        for(k=1;k<=NK[mglvl-1];k++)
      {
        p=i*2-1;
        q=j*2-1;
        r=k*2-1;
             for(itx=0;itx<=1;itx++)
             for(ity=0;ity<=1;ity++)
            for(itz=0;itz<=1;itz++)
            {
             IND=MIX(p+itx,q+ity,r+itz,mglvl);


               // var[IP][IND]= mginter(para, var[IP], fx[IND], fy[IND], fz[IND], i+itx-1, j+ity-1,  k+itz-1, mglvl-1);
                         //var[IP][IND]= var[IP][MIX(i,j,k,mglvl-1)];
              //var[IP][IND]=mginter_quad(para,var, var[IP], fx[IND], fy[IND], fz[IND], p+itx,q+ity,r+itz,i+itx-1, j+ity-1,  k+itz-1, mglvl-1);
            }
      }
}


void vel_map (PARA_DATA *para, REAL **var, int mglvl)
{

    int i,j,k;
  int p,q,r;
  int iter,it,jt,kt;
    int NMG=para->geom->NMG;
    int lvl;
    int NI[10],NJ[10],NK[10],IMAX[10],IJMAX[10],IJK[10];
  REAL *fx=var[FX],*fy=var[FY],*fz=var[FZ];
  REAL *flagp=var[FLAGP],*flagu=var[FLAGU],*flagv=var[FLAGV],*flagw=var[FLAGW];
  REAL ab[8];
    REAL DX,DY,DZ,ax,ay,az,ap;
    REAL *u=var[VX],*v=var[VY],*w=var[VZ],*ip=var[IP];
    REAL *x = var[X], *y = var[Y], *z = var[Z];
  REAL dt=para->mytime->dt;
  REAL massmax=0;
  REAL res=1;
  //REAL *pp=var[PP];
  REAL pp[8];

  for(lvl=1;lvl<=NMG;lvl++)
  {
   NI[lvl]=para->geom->NI[lvl];
   NJ[lvl]=para->geom->NJ[lvl]; 
   NK[lvl]=para->geom->NK[lvl]; 
   IMAX[lvl]=para->geom->IMAX[lvl];
   IJMAX[lvl]=para->geom->IJMAX[lvl];
   IJK[lvl]=para->geom->IJK[lvl];
  }

        /***********exterior velocity on fine grid correction********************/
  
  /*
  for(i=0;i<=NI[mglvl-1];i++)
  {
     p=i*2;
      for(q=0;q<=NJ[mglvl];q++)
        for(r=0;r<=NK[mglvl];r++)
      {

              if(flagu[MIX(p,q,r,mglvl)]>=0) u[MIX(p,q,r,mglvl)]=u[MIX(p,q,r,mglvl)];
              else u[MIX(p,q,r,mglvl)] = u[MIX(p,q,r,mglvl)]- dt*(ip[MIX(p+1,q,r,mglvl)]/(x[MIX(p+1,q,r,mglvl)]-x[MIX(p,q,r,mglvl)])-ip[MIX(p,q,r,mglvl)]/(x[MIX(p+1,q,r,mglvl)]-x[MIX(p,q,r,mglvl)]));
      }
  }

  for(j=0;j<=NJ[mglvl-1];j++)
  {
     q=j*2;
      for(p=0;p<=NI[mglvl];p++)
        for(r=0;r<=NK[mglvl];r++)
      {
             if (flagv[MIX(p,q,r,mglvl)]>=0) v[MIX(p,q,r,mglvl)] = v[MIX(p,q,r,mglvl)];
              else v[MIX(p,q,r,mglvl)] = v[MIX(p,q,r,mglvl)]- dt*(ip[MIX(p,q+1,r,mglvl)]/(y[MIX(p,q+1,r,mglvl)]-y[MIX(p,q,r,mglvl)])-ip[MIX(p,q,r,mglvl)]/(y[MIX(p,q+1,r,mglvl)]-y[MIX(p,q,r,mglvl)]));
      }
  }

  for(j=0;j<=NJ[mglvl-1];j++)
  {
     r=k*2;
      for(q=0;q<=NJ[mglvl];q++)
        for(p=0;p<=NI[mglvl];p++)
      {

             if (flagw[MIX(p,q,r,mglvl)]>=0) w[MIX(p,q,r,mglvl)] = w[MIX(p,q,r,mglvl)];
              else w[MIX(p,q,r,mglvl)] = w[MIX(p,q,r,mglvl)]- dt*(ip[MIX(p,q,r+1,mglvl)]/(z[MIX(p,q,r+1,mglvl)]-z[MIX(p,q,r,mglvl)])-ip[MIX(p,q,r,mglvl)]/(z[MIX(p,q,r+1,mglvl)]-z[MIX(p,q,r,mglvl)]));
      }
  }
  */


    for(i=1;i<=NI[mglvl-1];i++)
      for(j=1;j<=NJ[mglvl-1];j++)
          for(k=1;k<=NK[mglvl-1];k++)
        {
          p=i*2;
          q=j*2;
          r=k*2;

               for(jt=0;jt>=-1;jt--)
             for(kt=0;kt>=-1;kt--)
             {
               if(flagu[MIX(p,q+jt,r+kt,mglvl)]>=0) u[MIX(p,q+jt,r+kt,mglvl)]=u[MIX(p,q+jt,r+kt,mglvl)];
                         else u[MIX(p,q+jt,r+kt,mglvl)] = u[MIX(p,q+jt,r+kt,mglvl)] - dt*(ip[MIX(i+1,j,k,mglvl-1)]/(x[MIX(i+1,j,k,mglvl-1)]-x[MIX(i,j,k,mglvl-1)])-ip[MIX(i,j,k,mglvl-1)]/(x[MIX(i+1,j,k,mglvl-1)]-x[MIX(i,j,k,mglvl-1)]));
             
             }
             
             for(it=0;it>=-1;it--)
             for(kt=0;kt>=-1;kt--)
             {
               if(flagv[MIX(p+it,q,r+kt,mglvl)]>=0) v[MIX(p+it,q,r+kt,mglvl)]=v[MIX(p+it,q,r+kt,mglvl)];
                         else v[MIX(p+it,q,r+kt,mglvl)] = v[MIX(p+it,q,r+kt,mglvl)]- dt*(ip[MIX(i,j+1,k,mglvl-1)]/(y[MIX(i,j+1,k,mglvl-1)]-y[MIX(i,j,k,mglvl-1)])-ip[MIX(i,j,k,mglvl-1)]/(y[MIX(i,j+1,k,mglvl-1)]-y[MIX(i,j,k,mglvl-1)]));
             
             }

               for(jt=0;jt>=-1;jt--)
             for(it=0;it>=-1;it--)
             {
               if(flagw[MIX(p+it,q+jt,r,mglvl)]>=0) w[MIX(p+it,q+jt,r,mglvl)]=w[MIX(p+it,q+jt,r,mglvl)];
                         else w[MIX(p+it,q+jt,r,mglvl)] = w[MIX(p+it,q+jt,r,mglvl)]- dt*(ip[MIX(i,j,k+1,mglvl-1)]/(z[MIX(i,j,k+1,mglvl-1)]-z[MIX(i,j,k,mglvl-1)])-ip[MIX(i,j,k,mglvl-1)]/(z[MIX(i,j,k+1,mglvl-1)]-z[MIX(i,j,k,mglvl-1)]));
             
             }
        
            
        
        }


        /***********interior velocity on fine grid correction********************/
     
  for(i=1;i<=NI[mglvl-1];i++)
      for(j=1;j<=NJ[mglvl-1];j++)
        for(k=1;k<=NK[mglvl-1];k++)
      {

              if(flagp[MIX(i,j,k,mglvl-1)]>=0) continue;
        else
        {
        p=i*2-1;
        q=j*2-1;
        r=k*2-1;

          DX=(x[MIX(p+1,q,r,mglvl)]-x[MIX(p,q,r,mglvl)]);
          DY=(y[MIX(p,q+1,r,mglvl)]-y[MIX(p,q,r,mglvl)]);
          DZ=(z[MIX(p,q,r+1,mglvl)]-z[MIX(p,q,r,mglvl)]);

        ax=DY*DZ/DX;
        ay=DZ*DX/DY;
        az=DX*DY/DZ;
        ap=ax+ay+az;

        for(it=0;it<=1;it++)
          for(jt=0;jt<=1;jt++)
            for(kt=0;kt<=1;kt++)
            {

            ab[PIX(it,jt,kt)]=DX*DY*DZ*((u[MIX(p+it-1,q+jt,r+kt,mglvl)]-u[MIX(p+it,q+jt,r+kt,mglvl)])/DX
                                     +(v[MIX(p+it,q+jt-1,r+kt,mglvl)]-v[MIX(p+it,q+jt,r+kt,mglvl)])/DY
                         +(w[MIX(p+it,q+jt,r+kt-1,mglvl)]-w[MIX(p+it,q+jt,r+kt,mglvl)])/DZ);
            
            }

            /*
          c[GIX(0,0)]=ap;c[GIX(0,1)]=-ax;c[GIX(0,2)]=-ay;c[GIX(0,4)]=-az;
              c[GIX(1,1)]=ap;c[GIX(1,0)]=-ax;c[GIX(1,3)]=-ay;c[GIX(1,5)]=-az;
              c[GIX(2,2)]=ap;c[GIX(2,3)]=-ax;c[GIX(2,0)]=-ay;c[GIX(2,6)]=-az;
              c[GIX(3,3)]=ap;c[GIX(3,2)]=-ax;c[GIX(3,1)]=-ay;c[GIX(3,7)]=-az;
              c[GIX(4,4)]=ap;c[GIX(4,5)]=-ax;c[GIX(4,6)]=-ay;c[GIX(4,0)]=-az;
              c[GIX(5,5)]=ap;c[GIX(5,4)]=-ax;c[GIX(5,7)]=-ay;c[GIX(5,1)]=-az;
              c[GIX(6,6)]=ap;c[GIX(6,7)]=-ax;c[GIX(6,4)]=-ay;c[GIX(6,2)]=-az;
              c[GIX(7,7)]=1.0f;
        ab[7]=0;

         direct_solver( c, pp, ab);
         */

            
               for(it=0;it<=7;it++) pp[it]=0;



      /*  pp[MIX(p,q,r,mglvl)]=(ax*pp[MIX(p+1,q,r,mglvl)]+ay*pp[MIX(p,q+1,r,mglvl)]+az*pp[MIX(p,q,r+1,mglvl)]+ab[0])/ap;
              pp[MIX(p+1,q,r,mglvl)]=(ax*pp[MIX(p,q,r,mglvl)]+ay*pp[MIX(p+1,q+1,r,mglvl)]+az*pp[MIX(p+1,q,r+1,mglvl)]+ab[1])/ap;
              pp[MIX(p,q+1,r,mglvl)]=(ax*pp[MIX(p+1,q+1,r,mglvl)]+ay*pp[MIX(p,q,r,mglvl)]+az*pp[MIX(p,q+1,r+1,mglvl)]+ab[2])/ap;
              pp[MIX(p+1,q+1,r,mglvl)]=(ax*pp[MIX(p,q+1,r,mglvl)]+ay*pp[MIX(p+1,q,r,mglvl)]+az*pp[MIX(p+1,q+1,r+1,mglvl)]+ab[3])/ap;
              pp[MIX(p,q,r+1,mglvl)]=(ax*pp[MIX(p+1,q,r+1,mglvl)]+ay*pp[MIX(p,q+1,r+1,mglvl)]+az*pp[MIX(p,q,r,mglvl)]+ab[4])/ap;
              pp[MIX(p+1,q,r+1,mglvl)]=(ax*pp[MIX(p,q,r+1,mglvl)]+ay*pp[MIX(p+1,q+1,r+1,mglvl)]+az*pp[MIX(p+1,q,r,mglvl)]+ab[5])/ap;
              pp[MIX(p,q+1,r+1,mglvl)]=(ax*pp[MIX(p+1,q+1,r+1,mglvl)]+ay*pp[MIX(p,q,r+1,mglvl)]+az*pp[MIX(p,q+1,r,mglvl)]+ab[6])/ap;
              pp[MIX(p+1,q+1,r+1,mglvl)]=(ax*pp[MIX(p,q+1,r+1,mglvl)]+ay*pp[MIX(p+1,q,r+1,mglvl)]+az*pp[MIX(p+1,q+1,r,mglvl)]+ab[7])/ap;
        }

                      it=0;
        for(jt=0;jt<=1;jt++)
            for(kt=0;kt<=1;kt++)
              u[MIX(p+it,q+jt,r+kt,mglvl)] -= (pp[MIX(p+1,q+jt,r+kt,mglvl)]-pp[MIX(p,q+jt,r+kt,mglvl)])/DX;

              jt=0;
        for(it=0;it<=1;it++)
            for(kt=0;kt<=1;kt++)
              v[MIX(p+it,q+jt,r+kt,mglvl)] -= (pp[MIX(p+it,q+1,r+kt,mglvl)]-pp[MIX(p+it,q,r+kt,mglvl)])/DY;

              kt=0;
        for(jt=0;jt<=1;jt++)
            for(it=0;it<=1;it++)
              w[MIX(p+it,q+jt,r+kt,mglvl)] -= (pp[MIX(p+it,q+jt,r+1,mglvl)]-pp[MIX(p+it,q+jt,r,mglvl)])/DZ;

           */
          iter=0;
              while(iter<10) //&& res>0.000001
        {
          iter++;
          res=0;

         pp[0]=(ax*pp[1]+ay*pp[2]+az*pp[4]+ab[0])/ap;
              pp[1]=(ax*pp[0]+ay*pp[3]+az*pp[5]+ab[1])/ap;
              pp[2]=(ax*pp[3]+ay*pp[0]+az*pp[6]+ab[2])/ap;
              pp[3]=(ax*pp[2]+ay*pp[1]+az*pp[7]+ab[3])/ap;
              pp[4]=(ax*pp[5]+ay*pp[6]+az*pp[0]+ab[4])/ap;
              pp[5]=(ax*pp[4]+ay*pp[7]+az*pp[1]+ab[5])/ap;
        pp[6]=(ax*pp[7]+ay*pp[4]+az*pp[2]+ab[6])/ap;
              pp[7]=(ax*pp[6]+ay*pp[5]+az*pp[3]+ab[7])/ap;
        pp[6]=(ax*pp[7]+ay*pp[4]+az*pp[2]+ab[6])/ap;
              pp[5]=(ax*pp[4]+ay*pp[7]+az*pp[1]+ab[5])/ap;
              pp[4]=(ax*pp[5]+ay*pp[6]+az*pp[0]+ab[4])/ap;
              pp[3]=(ax*pp[2]+ay*pp[1]+az*pp[7]+ab[3])/ap;
              pp[2]=(ax*pp[3]+ay*pp[0]+az*pp[6]+ab[2])/ap;
              pp[1]=(ax*pp[0]+ay*pp[3]+az*pp[5]+ab[1])/ap;
        pp[0]=(ax*pp[1]+ay*pp[2]+az*pp[4]+ab[0])/ap;

      /*  resp[0]=fabs(ap*pp[0]-(ax*pp[1]+ay*pp[2]+az*pp[4]+ab[0]));
       resp[1]=fabs(ap*pp[1]-(ax*pp[0]+ay*pp[3]+az*pp[5]+ab[1]));
        resp[2]=fabs(ap*pp[2]-(ax*pp[3]+ay*pp[0]+az*pp[6]+ab[2]));
       resp[3]=fabs(ap*pp[3]-(ax*pp[2]+ay*pp[1]+az*pp[7]+ab[3]));
        resp[4]=fabs(ap*pp[4]-(ax*pp[5]+ay*pp[6]+az*pp[0]+ab[4]));
       resp[5]=fabs(ap*pp[5]-(ax*pp[4]+ay*pp[7]+az*pp[1]+ab[5]));
       resp[6]=fabs(ap*pp[6]-(ax*pp[7]+ay*pp[4]+az*pp[2]+ab[6]));
       resp[7]=fabs(ap*pp[7]-(ax*pp[5]+ay*pp[4]+az*pp[2]+ab[6]));
              for(it=0;it<=6;it++) res=max(res,resp[it]);
        */
        }

                            it=0;
        for(jt=0;jt<=1;jt++)
            for(kt=0;kt<=1;kt++)
              u[MIX(p+it,q+jt,r+kt,mglvl)] -= (pp[PIX(1,jt,kt)]-pp[PIX(0,jt,kt)])/DX;

              jt=0;
        for(it=0;it<=1;it++)
            for(kt=0;kt<=1;kt++)
              v[MIX(p+it,q+jt,r+kt,mglvl)] -= (pp[PIX(it,1,kt)]-pp[PIX(it,0,kt)])/DY;

              kt=0;
        for(jt=0;jt<=1;jt++)
            for(it=0;it<=1;it++)
              w[MIX(p+it,q+jt,r+kt,mglvl)] -= (pp[PIX(it,jt,1)]-pp[PIX(it,jt,0)])/DZ;

         

    /*   
        

      resp[0]=fabs(ap*pp[0]-(ax*pp[1]+ay*pp[2]+az*pp[4]+ab[0]));
       resp[1]=fabs(ap*pp[1]-(ax*pp[0]+ay*pp[3]+az*pp[5]+ab[1]));
        resp[2]=fabs(ap*pp[2]-(ax*pp[3]+ay*pp[0]+az*pp[6]+ab[2]));
       resp[3]=fabs(ap*pp[3]-(ax*pp[2]+ay*pp[1]+az*pp[7]+ab[3]));
        resp[4]=fabs(ap*pp[4]-(ax*pp[5]+ay*pp[6]+az*pp[0]+ab[4]));
       resp[5]=fabs(ap*pp[5]-(ax*pp[4]+ay*pp[7]+az*pp[1]+ab[5]));
       resp[6]=fabs(ap*pp[6]-(ax*pp[7]+ay*pp[4]+az*pp[2]+ab[6]));
              for(iter=0;iter<=6;iter++) res=max(res,resp[iter]);
        */
    
        


        /*
        u[MIX(p,q,r,mglvl)]     -= (pp[1]-pp[0])/DX;
        u[MIX(p,q+1,r,mglvl)]   -= (pp[3]-pp[2])/DX;
        u[MIX(p,q,r+1,mglvl)]   -= (pp[5]-pp[4])/DX;
        u[MIX(p,q+1,r+1,mglvl)] -= (pp[7]-pp[6])/DX;

        v[MIX(p,q,r,mglvl)] -= (pp[2]-pp[0])/DY;
        v[MIX(p+1,q,r,mglvl)] -= (pp[3]-pp[1])/DY;
        v[MIX(p,q,r+1,mglvl)] -= (pp[6]-pp[4])/DY;
        v[MIX(p+1,q,r+1,mglvl)] -= (pp[7]-pp[5])/DY;

        w[MIX(p,q,r,mglvl)] -= (pp[4]-pp[0])/DZ;
        w[MIX(p+1,q,r,mglvl)] -= (pp[5]-pp[1])/DZ;
        w[MIX(p,q+1,r,mglvl)] -= (pp[6]-pp[2])/DZ;
        w[MIX(p+1,q+1,r,mglvl)] -= (pp[7]-pp[3])/DZ;
        */


             // for(it=0;it<=1;it++)
        //  for(jt=0;jt<=1;jt++)
        //    for(kt=0;kt<=1;kt++)
        //    {
            //     resmass= fabs((u[MIX(p+it-1,q+jt,r+kt,mglvl)]-u[MIX(p+it,q+jt,r+kt,mglvl)])*DY*DZ
                  //        + (v[MIX(p+it,q+jt-1,r+kt,mglvl)]-v[MIX(p+it,q+jt,r+kt,mglvl)])*DX*DZ
        //      + (w[MIX(p+it-1,q+jt,r+kt-1,mglvl)]-w[MIX(p+it,q+jt,r+kt,mglvl)])*DY*DX);

          //   if(massmax<resmass) massmax=resmass;
          
          //  }
          
                     // printf("%f\t",massmax);
          }
               

      } //end_loop
      
    
}

void cmap (PARA_DATA *para, REAL **var, int mglvl)
{
    int i,j,k;
  int p,q,r;
    int NMG=para->geom->NMG;
  REAL *u = var[VX], *v = var[VY], *w = var[VZ];
    int lvl;
    int NI[10],NJ[10],NK[10],IMAX[10],IJMAX[10],IJK[10];

  for(lvl=1;lvl<=NMG;lvl++)
  {
   NI[lvl]=para->geom->NI[lvl];
   NJ[lvl]=para->geom->NJ[lvl]; 
   NK[lvl]=para->geom->NK[lvl]; 
   IMAX[lvl]=para->geom->IMAX[lvl];
   IJMAX[lvl]=para->geom->IJMAX[lvl];
   IJK[lvl]=para->geom->IJK[lvl];

  }


  for(i=0;i<=NI[mglvl];i++)
      for(j=1;j<=NJ[mglvl];j++)
        for(k=1;k<=NK[mglvl];k++)
      {
        p=i*2;
        q=j*2;
        r=k*2;
          u[MIX(i,j,k,mglvl)]= (u[MIX(p,q,r,mglvl+1)]  +u[MIX(p,q,r-1,mglvl+1)]  +u[MIX(p,q-1,r,mglvl+1)]  +u[MIX(p,q-1,r-1,mglvl+1)])/4.0f;
      }

  for(i=1;i<=NI[mglvl];i++)
      for(j=0;j<=NJ[mglvl];j++)
        for(k=1;k<=NK[mglvl];k++)
      {
        p=i*2;
        q=j*2;
        r=k*2;
          v[MIX(i,j,k,mglvl)]= (v[MIX(p,q,r,mglvl+1)]  +v[MIX(p,q,r-1,mglvl+1)]  +v[MIX(p-1,q,r,mglvl+1)]  +v[MIX(p-1,q,r-1,mglvl+1)])/4.0f;
      }

  for(i=1;i<=NI[mglvl];i++)
      for(j=1;j<=NJ[mglvl];j++)
        for(k=0;k<=NK[mglvl];k++)
      {
        p=i*2;
        q=j*2;
        r=k*2;
          w[MIX(i,j,k,mglvl)]= (w[MIX(p,q,r,mglvl+1)]  +w[MIX(p,q-1,r,mglvl+1)]  +w[MIX(p-1,q,r,mglvl+1)]  +w[MIX(p-1,q-1,r,mglvl+1)])/4.0f;
      }


}