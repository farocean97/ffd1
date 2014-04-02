#include <math.h>
#include <stdio.h>
#include "data_structure.h"
#include "boundary.h"
#include "advection.h"
#include "solver.h"
#include "utility.h"
#include "interpolation.h"

/******************************************************************************
| advection
******************************************************************************/
void traceback(PARA_DATA *para, REAL **var, int **BINDEX){
  int i, j, k;
  int imax = para->geom->imax, jmax = para->geom->jmax;
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  REAL *fx=var[FX],*fy=var[FY],*fz=var[FZ];
  REAL dt = para->mytime->dt; 
  REAL dx,dy,dz;
  REAL u0, v0, w0;
  REAL *x = var[X],  *y = var[Y],  *z = var[Z]; 
  REAL *gx = var[GX],  *gy = var[GY],  *gz = var[GZ]; 
  REAL *u = var[VX], *v = var[VY], *w = var[VZ];
  REAL *flagp = var[FLAGP],*flagu = var[FLAGU];
  REAL *flagv = var[FLAGV],*flagw = var[FLAGW];
  int  COOD[3], LOC[3];
  REAL OL[3];
  int  OC[3];
  int count;

  count=para->mytime->t_step%6;

  FOR_EACH_CELL
    if(flagp[FIX(i,j,k)]>0) continue;
    u0 = 0.5f*(u[FIX(i,j,k)]+u[FIX(i-1,j,k)]);
    v0 = 0.5f*(v[FIX(i,j,k)]+v[FIX(i,j-1,k)]);
    w0 = 0.5f*(w[FIX(i,j,k)]+w[FIX(i,j,k-1)]);
          
    OL[X] = x[FIX(i,j,k)] - u0*dt; 
    OL[Y] = y[FIX(i,j,k)] - v0*dt;
    OL[Z] = z[FIX(i,j,k)] - w0*dt;
    
    OC[X] = i; OC[Y] = j; OC[Z] = k;  
    COOD[X] =1; COOD[Y]=1; COOD[Z]=1;
    LOC[X]  =1; LOC[Y] =1; LOC[Z] =1;

    while(COOD[X]==1 || COOD[Y] ==1 || COOD[Z] == 1){
     switch(count) {
        case 0:
          if(COOD[X]==1 && LOC[X]==1)
            XLOCATION(para, var, flagp, x, u0, i, j, k, OL,OC, LOC ,COOD);
          if(COOD[Y]==1 && LOC[Y]==1)
            YLOCATION(para, var, flagp, y, v0, i, j, k, OL,OC, LOC ,COOD); 
          if(COOD[Z]==1 && LOC[Z]==1)
            ZLOCATION(para, var, flagp, z, w0, i, j, k, OL,OC, LOC ,COOD); 
          break;
        case 1:
          if(COOD[X]==1 && LOC[X]==1)
            XLOCATION(para, var, flagp, x, u0, i, j, k, OL,OC, LOC ,COOD);
          if(COOD[Z]==1 && LOC[Z]==1)
            ZLOCATION(para, var, flagp, z, w0, i, j, k, OL,OC, LOC ,COOD); 
          if(COOD[Y]==1 && LOC[Y]==1)
            YLOCATION(para, var, flagp, y, v0, i, j, k, OL,OC, LOC ,COOD); 
          break;
        case 2:
          if(COOD[Y]==1 && LOC[Y]==1)
            YLOCATION(para, var, flagp, y, v0, i, j, k, OL,OC, LOC ,COOD); 
          if(COOD[X]==1 && LOC[X]==1)
            XLOCATION(para, var, flagp, x, u0, i, j, k, OL,OC, LOC ,COOD);
          if(COOD[Z]==1 && LOC[Z]==1)
            ZLOCATION(para, var, flagp, z, w0, i, j, k, OL,OC, LOC ,COOD); 
          break;
        case 3:
          if(COOD[Y]==1 && LOC[Y]==1)
            YLOCATION(para, var, flagp, y, v0, i, j, k, OL,OC, LOC ,COOD); 
          if(COOD[Z]==1 && LOC[Z]==1)
            ZLOCATION(para, var, flagp, z, w0, i, j, k, OL,OC, LOC ,COOD); 
          if(COOD[X]==1 && LOC[X]==1)
            XLOCATION(para, var, flagp, x, u0, i, j, k, OL,OC, LOC ,COOD);
          break;
        case 4:
          if(COOD[Z]==1 && LOC[Z]==1)
            ZLOCATION(para, var, flagp, z, w0, i, j, k, OL,OC, LOC ,COOD); 
          if(COOD[X]==1 && LOC[X]==1)
            XLOCATION(para, var, flagp, x, u0, i, j, k, OL,OC, LOC ,COOD);
          if(COOD[Y]==1 && LOC[Y]==1)
            YLOCATION(para, var, flagp, y, v0, i, j, k, OL,OC, LOC ,COOD); 
          break;
        case 5:
          if(COOD[Z]==1 && LOC[Z]==1)
            ZLOCATION(para, var, flagp, z, w0, i, j, k, OL,OC, LOC ,COOD); 
          if(COOD[Y]==1 && LOC[Y]==1)
            YLOCATION(para, var, flagp, y, v0, i, j, k, OL,OC, LOC ,COOD); 
          if(COOD[X]==1 && LOC[X]==1)
            XLOCATION(para, var, flagp, x, u0, i, j, k, OL,OC, LOC ,COOD);
          break;
     }
      
    } 
 
     if(u0>=0 && LOC[X] == 0) OC[X] -=1;
	   if(v0>=0 && LOC[Y] == 0) OC[Y] -=1;
	   if(w0>=0 && LOC[Z] == 0) OC[Z] -=1;

    if(u0<0 && LOC[X]>0) OC[X] -=1;
    if(v0<0 && LOC[Y]>0) OC[Y] -=1;
    if(w0<0 && LOC[Z]>0) OC[Z] -=1;

   
    BINDEX[4][FIX(i,j,k)]=OC[X];
    BINDEX[5][FIX(i,j,k)]=OC[Y];
    BINDEX[6][FIX(i,j,k)]=OC[Z];


    fx[FIX(i,j,k)] = (OL[X]- x[FIX(OC[X],OC[Y],OC[Z])])
                    /(x[FIX(OC[X]+1,OC[Y],OC[Z])]- x[FIX(OC[X],OC[Y],OC[Z])]); 
    fy[FIX(i,j,k)] = (OL[Y]- y[FIX(OC[X],OC[Y],OC[Z])])
                    /(y[FIX(OC[X],OC[Y]+1,OC[Z])]- y[FIX(OC[X],OC[Y],OC[Z])]);
    fz[FIX(i,j,k)] = (OL[Z]- z[FIX(OC[X],OC[Y],OC[Z])])
                    /(z[FIX(OC[X],OC[Y],OC[Z]+1)]- z[FIX(OC[X],OC[Y],OC[Z])]);

  
  END_FOR

}

void traceback_UVW(PARA_DATA *para, REAL **var, int var_type, REAL *flag, int **BINDEX){
  int i, j, k;
  int imax = para->geom->imax, jmax = para->geom->jmax;
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  REAL *fx=var[FX],*fy=var[FY],*fz=var[FZ];
  REAL dt = para->mytime->dt; 
  REAL dx,dy,dz;
  REAL u0, v0, w0;
  REAL *x = var[X],  *y = var[Y],  *z = var[Z]; 
  REAL *gx = var[GX],  *gy = var[GY],  *gz = var[GZ]; 
  REAL *u = var[VX], *v = var[VY], *w = var[VZ];
  int  COOD[3], LOC[3];
  REAL OL[3];
  int  OC[3];

  switch(var_type) {
    case VX:
      FOR_EACH_CELL
        if(flag[FIX(i,j,k)]>=0) continue;
        u0 = u[FIX(i,j,k)];
        v0 = 0.5f *( (v[FIX(i,j,k)]+ v[FIX(i,j-1,k)])
                    *(x[FIX(i+1,j,k)]-gx[FIX(i,j,k)])
		                +(v[FIX(i+1,j,k)]+v[FIX(i+1,j-1,k)])
                    *(gx[FIX(i,  j,k)]-x[FIX(i,j,k)])) 
                    /(x[FIX(i+1,j,k)]-x[FIX(i,j,k)]);
	
	      w0 = 0.5f *( (w[FIX(i,j,k)]+ w[FIX(i,j,k-1)])
                    *(x[FIX(i+1,j,k)]-gx[FIX(i,j,k)])
		                +(w[FIX(i+1,j,k)]+ w[FIX(i+1,j, k-1)])
                    *(gx[FIX(i,j,k)]-x[FIX(i,j,k)]))
                    /(x[FIX(i+1,j,k)]-x[FIX(i,j,k)]); 
          
        OL[X] = gx[FIX(i,j,k)] - u0*dt; 
        OL[Y] = y[FIX(i,j,k)] - v0*dt;
        OL[Z] = z[FIX(i,j,k)] - w0*dt;
    
        OC[X] = i; OC[Y] = j; OC[Z] = k;  
        COOD[X] =1; COOD[Y]=1; COOD[Z]=1;
        LOC[X]  =1; LOC[Y] =1; LOC[Z] =1;

        while(COOD[X]==1 || COOD[Y] ==1 || COOD[Z] == 1){
          if(COOD[X]==1 && LOC[X]==1)
            XLOCATION_U(para, var, flag, gx, u0, i, j, k, OL,OC, LOC ,COOD);
          if(COOD[Y]==1 && LOC[Y]==1)
            YLOCATION(para, var, flag, y, v0, i, j, k, OL,OC, LOC ,COOD); 
          if(COOD[Z]==1 && LOC[Z]==1)
            ZLOCATION(para, var, flag, z, w0, i, j, k, OL,OC, LOC ,COOD); 
        }
 	      if(u0>=0 && LOC[X] == 0) OC[X] -=1;
	      if(v0>=0 && LOC[Y] == 0) OC[Y] -=1;
	      if(w0>=0 && LOC[Z] == 0) OC[Z] -=1;

        if(u0<0 && LOC[X]>0) OC[X] -=1;
        if(v0<0 && LOC[Y]>0) OC[Y] -=1;
        if(w0<0 && LOC[Z]>0) OC[Z] -=1;

   
        BINDEX[4][FIX(i,j,k)]=OC[X];
        BINDEX[5][FIX(i,j,k)]=OC[Y];
        BINDEX[6][FIX(i,j,k)]=OC[Z];

        fx[FIX(i,j,k)] = (OL[X]- gx[FIX(OC[X],OC[Y],OC[Z])])
                    /(gx[FIX(OC[X]+1,OC[Y],OC[Z])]- gx[FIX(OC[X],OC[Y],OC[Z])]); 
        fy[FIX(i,j,k)] = (OL[Y]- y[FIX(OC[X],OC[Y],OC[Z])])
                    /(y[FIX(OC[X],OC[Y]+1,OC[Z])]- y[FIX(OC[X],OC[Y],OC[Z])]);
        fz[FIX(i,j,k)] = (OL[Z]- z[FIX(OC[X],OC[Y],OC[Z])])
                    /(z[FIX(OC[X],OC[Y],OC[Z]+1)]- z[FIX(OC[X],OC[Y],OC[Z])]); 
      END_FOR
      break;
    case VY:
      FOR_EACH_CELL
        if(flag[FIX(i,j,k)]>=0) continue;
	        u0 = 0.5f* ( (u[FIX(i,j,k)]+u[FIX(i-1,j,k)])
                      *(y[FIX(i,j+1,k)]-gy[FIX(i,j,k)])
                      +(u[FIX(i,j+1,k)]+u[FIX(i-1,j+1,k)])
                      *(gy[FIX(i,j,k)]-y[FIX(i,j,k)]))
                      /(y[FIX(i,j+1,k)]-y[FIX(i,j,k)]);
	
          v0 = v[FIX(i,j,k)]; 

        	w0 = 0.5f *( (w[FIX(i,j,k)]+ w[FIX(i,j,k-1)])
                      *(y[FIX(i,j+1,k)]- gy[FIX(i,j,k)])
                      +(w[FIX(i,j+1,k)]+ w[FIX(i,j+1,k-1)])
                      *(gy[FIX(i,j,k)]-y[FIX(i,j,k)]))
                      /(y[FIX(i,j+1,k)]-y[FIX(i,j,k)]); 
          
        OL[X] = x[FIX(i,j,k)] - u0*dt; 
        OL[Y] = gy[FIX(i,j,k)] - v0*dt;
        OL[Z] = z[FIX(i,j,k)] - w0*dt;
    
        OC[X] = i; OC[Y] = j; OC[Z] = k;  
        COOD[X] =1; COOD[Y]=1; COOD[Z]=1;
        LOC[X]  =1; LOC[Y] =1; LOC[Z] =1;

        while(COOD[X]==1 || COOD[Y] ==1 || COOD[Z] == 1){
          if(COOD[X]==1 && LOC[X]==1)
            XLOCATION(para, var, flag, x, u0, i, j, k, OL,OC, LOC ,COOD);
          if(COOD[Y]==1 && LOC[Y]==1)
            YLOCATION_V(para, var, flag, gy, v0, i, j, k, OL,OC, LOC ,COOD); 
          if(COOD[Z]==1 && LOC[Z]==1)
            ZLOCATION(para, var, flag, z, w0, i, j, k, OL,OC, LOC ,COOD); 
        }
 	      if(u0>=0 && LOC[X] == 0) OC[X] -=1;
	      if(v0>=0 && LOC[Y] == 0) OC[Y] -=1;
	      if(w0>=0 && LOC[Z] == 0) OC[Z] -=1;

        if(u0<0 && LOC[X]>0) OC[X] -=1;
        if(v0<0 && LOC[Y]>0) OC[Y] -=1;
        if(w0<0 && LOC[Z]>0) OC[Z] -=1;

   
        BINDEX[4][FIX(i,j,k)]=OC[X];
        BINDEX[5][FIX(i,j,k)]=OC[Y];
        BINDEX[6][FIX(i,j,k)]=OC[Z];

        fx[FIX(i,j,k)] = (OL[X]- x[FIX(OC[X],OC[Y],OC[Z])])
                    /(x[FIX(OC[X]+1,OC[Y],OC[Z])]- x[FIX(OC[X],OC[Y],OC[Z])]); 
        fy[FIX(i,j,k)] = (OL[Y]- gy[FIX(OC[X],OC[Y],OC[Z])])
                    /(gy[FIX(OC[X],OC[Y]+1,OC[Z])]- gy[FIX(OC[X],OC[Y],OC[Z])]);
        fz[FIX(i,j,k)] = (OL[Z]- z[FIX(OC[X],OC[Y],OC[Z])])
                    /(z[FIX(OC[X],OC[Y],OC[Z]+1)]- z[FIX(OC[X],OC[Y],OC[Z])]); 
      END_FOR
      break;
    case VZ:
      FOR_EACH_CELL
        if(flag[FIX(i,j,k)]>=0) continue;
	        u0 = 0.5f*( (u[FIX(i,j,k)]+ u[FIX(i-1,j,k)])
                     *(z[FIX(i,j,k+1)]-gz[FIX(i,j,k)])
                     +(u[FIX(i,j,k+1)]+ u[FIX(i-1,j,k+1)])
                     *(gz[FIX(i,j,k )]- z[FIX(i,j,k)]))
                     /(z[FIX(i,j,k+1)]-z[FIX(i,j,k)]);

	        v0 = 0.5f*( (v[FIX(i,j,k)]+ v[FIX(i,j-1,k)])
                     *(z[FIX(i,j,k+1)]-gz[FIX(i,j,k)])
                     +(v[FIX(i,j,k+1)]+v[FIX(i,j-1,k+1)])
                     *(gz[FIX(i,j,k)]-z [FIX(i,j,k)]))
                     /(z[FIX(i,j,k+1)]-z[FIX(i,j,k)]); 

	        w0 = w[FIX(i,j,k)]; 
          
          OL[X] = x[FIX(i,j,k)] - u0*dt; 
          OL[Y] = y[FIX(i,j,k)] - v0*dt;
          OL[Z] = gz[FIX(i,j,k)] - w0*dt;
    
          OC[X] = i; OC[Y] = j; OC[Z] = k;  
          COOD[X] =1; COOD[Y]=1; COOD[Z]=1;
          LOC[X]  =1; LOC[Y] =1; LOC[Z] =1;

          while(COOD[X]==1 || COOD[Y] ==1 || COOD[Z] == 1){
            if(COOD[X]==1 && LOC[X]==1)
              XLOCATION(para, var, flag, x, u0, i, j, k, OL,OC, LOC ,COOD);
            if(COOD[Y]==1 && LOC[Y]==1)
              YLOCATION(para, var, flag, y, v0, i, j, k, OL,OC, LOC ,COOD); 
            if(COOD[Z]==1 && LOC[Z]==1)
              ZLOCATION_W(para, var, flag, gz, w0, i, j, k, OL,OC, LOC ,COOD); 
          }
 	        if(u0>=0 && LOC[X] == 0) OC[X] -=1;
	        if(v0>=0 && LOC[Y] == 0) OC[Y] -=1;
	        if(w0>=0 && LOC[Z] == 0) OC[Z] -=1;

          if(u0<0 && LOC[X]>0) OC[X] -=1;
          if(v0<0 && LOC[Y]>0) OC[Y] -=1;
          if(w0<0 && LOC[Z]>0) OC[Z] -=1;

   
          BINDEX[4][FIX(i,j,k)]=OC[X];
          BINDEX[5][FIX(i,j,k)]=OC[Y];
          BINDEX[6][FIX(i,j,k)]=OC[Z];

          fx[FIX(i,j,k)] = (OL[X]- x[FIX(OC[X],OC[Y],OC[Z])])
                    /(x[FIX(OC[X]+1,OC[Y],OC[Z])]- x[FIX(OC[X],OC[Y],OC[Z])]); 
          fy[FIX(i,j,k)] = (OL[Y]- y[FIX(OC[X],OC[Y],OC[Z])])
                    /(y[FIX(OC[X],OC[Y]+1,OC[Z])]- y[FIX(OC[X],OC[Y],OC[Z])]);
          fz[FIX(i,j,k)] = (OL[Z]- gz[FIX(OC[X],OC[Y],OC[Z])])
                    /(gz[FIX(OC[X],OC[Y],OC[Z]+1)]- gz[FIX(OC[X],OC[Y],OC[Z])]); 
        END_FOR
        break;
      }
} 


void advection(PARA_DATA *para, REAL **var, int var_type,
                REAL *flag, REAL *d, REAL *d0, int **BINDEX){
  int i, j, k;
  int imax = para->geom->imax, jmax = para->geom->jmax;
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  int p,q,r;
  REAL *x = var[X], *y = var[Y], *z = var[Z];
  REAL *u=var[VX],*v=var[VY],*w=var[VZ];
  REAL *gx = var[GX], *gy = var[GY], *gz = var[GZ];
  REAL *vis=var[DEN];
  REAL dt = para->mytime->dt; 
  REAL *fx=var[FX],*fy=var[FY],*fz=var[FZ];

  switch(var_type){
    case TEMP:
    case DEN:
      traceback(para,var,BINDEX);
      FOR_EACH_CELL
        if(flag[FIX(i,j,k)]==1) continue;
        else if(flag[FIX(i,j,k)]==0){
          if(para->prob->plume_mod==1) {
            //d[FIX(i,j,k)]= d0[FIX(i,j,k)];
            //d[FIX(i,j,k)]= advection_upwind(para,var,d0,i,j,k);

            p=BINDEX[4][FIX(i,j,k)];
            q=BINDEX[5][FIX(i,j,k)];
            r=BINDEX[6][FIX(i,j,k)];
            d[FIX(i,j,k)] = interpolation_temp(para,var, d0, fx[FIX(i,j,k)], 
                                  fy[FIX(i,j,k)], fz[FIX(i,j,k)],i,j,k, p,q,r);
          }
          else {
            d[FIX(i,j,k)]= advection_upwind(para,var,d0,i,j,k);
          }

        }
        else
        {
          p=BINDEX[4][FIX(i,j,k)];
          q=BINDEX[5][FIX(i,j,k)];
          r=BINDEX[6][FIX(i,j,k)];

          d[FIX(i,j,k)] = interpolation_temp(para,var, d0, fx[FIX(i,j,k)], 
                                  fy[FIX(i,j,k)], fz[FIX(i,j,k)],i,j,k, p,q,r);
          //if(i==12 && j==9 && k==10)

 /*         if(i==9 && j==5 && k==8) printf("pqr=%d\t%d\t%d\t%f\t%f\t%f\t%f\n",p,q,r, fx[FIX(i,j,k)], fy[FIX(i,j,k)], fz[FIX(i,j,k)],d[FIX(i,j,k)]);
          if(i==9 && j==5 && k==8) printf("pqr=%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",d0[FIX(p,q,r)], d0[FIX(p,q,r+1)], d0[FIX(p,q+1,r)],d0[FIX(p,q+1,r+1)],
                                                                                        d0[FIX(p+1,q,r)], d0[FIX(p+1,q,r+1)], d0[FIX(p+1,q+1,r)],d0[FIX(p+1,q+1,r+1)]);*/
        }
      END_FOR
      break;
    case VX:
      traceback_UVW(para, var, var_type, flag, BINDEX);
      FOR_EACH_CELL
        if(flag[FIX(i,j,k)]>=0) continue;
        p=BINDEX[4][FIX(i,j,k)];
        q=BINDEX[5][FIX(i,j,k)];
        r=BINDEX[6][FIX(i,j,k)];
        d[FIX(i,j,k)] = interpolation(para, d0, fx[FIX(i,j,k)], 
                                 fy[FIX(i,j,k)], fz[FIX(i,j,k)], p,q,r);
      END_FOR
      break;
    case VY:
      traceback_UVW(para, var, var_type, flag, BINDEX);
      FOR_EACH_CELL
        if(flag[FIX(i,j,k)]>=0) continue;
        p=BINDEX[4][FIX(i,j,k)];
        q=BINDEX[5][FIX(i,j,k)];
        r=BINDEX[6][FIX(i,j,k)];
        d[FIX(i,j,k)] = interpolation(para, d0, fx[FIX(i,j,k)],
                                  fy[FIX(i,j,k)], fz[FIX(i,j,k)], p,q,r);
      END_FOR
      break;
    case VZ:
      traceback_UVW(para, var, var_type, flag, BINDEX);
      FOR_EACH_CELL
        if(flag[FIX(i,j,k)]>=0) continue;
        p=BINDEX[4][FIX(i,j,k)];
        q=BINDEX[5][FIX(i,j,k)];
        r=BINDEX[6][FIX(i,j,k)];
        d[FIX(i,j,k)] = interpolation(para, d0, fx[FIX(i,j,k)], 
                                  fy[FIX(i,j,k)], fz[FIX(i,j,k)], p,q,r);
      END_FOR
      break;
   }
      
}
  
void XLOCATION(PARA_DATA *para, REAL **var, REAL *flag, REAL *x, REAL u0,
               int i, int j, int k,  REAL *OL, int *OC, int *LOC , int *COOD){

  int imax = para->geom->imax, jmax = para->geom->jmax; 
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  REAL *u=var[VX];
  REAL dt=para->mytime->dt;
  REAL *flagu=var[FLAGU];

  if(OL[X]==x[FIX(OC[X],OC[Y],OC[Z])]){
    COOD[X]=0;
  }
  else if(OL[X]<x[FIX(OC[X],OC[Y],OC[Z])])
  {
    if(OC[X]>0) OC[X] -=1;
    if(OL[X]>=x[FIX(OC[X],OC[Y],OC[Z])]) COOD[X]=0;
    if(flag[FIX(OC[X],OC[Y],OC[Z])]>0 || flagu[FIX(OC[X],OC[Y],OC[Z])]==4){
      if(flagu[FIX(OC[X],OC[Y],OC[Z])]==0){
       if(OL[X]<x[FIX(OC[X],OC[Y],OC[Z])]) OL[X]=x[FIX(OC[X],OC[Y],OC[Z])];
      }
      else {
        OL[X]=x[FIX(OC[X]+1,OC[Y],OC[Z])];
      }
      LOC[X]=0;
      COOD[X]=0;
      OC[X] +=1;
    }
  }
  else {
    if(OC[X]<=imax) OC[X] +=1;
    if(OL[X] <=x[FIX(OC[X],OC[Y],OC[Z])]) COOD[X]=0;
    if(flag[FIX(OC[X],OC[Y],OC[Z])]>0 || flagu[FIX(OC[X]-1,OC[Y],OC[Z])]==4){
      if(flagu[FIX(OC[X]-1,OC[Y],OC[Z])]==0){
        if(OL[X]>x[FIX(OC[X],OC[Y],OC[Z])]) OL[X]=x[FIX(OC[X],OC[Y],OC[Z])];
      }
      else {
        OL[X]=x[FIX(OC[X]-1,OC[Y],OC[Z])];
      }
      LOC[X]=0;
      COOD[X]=0;
      OC[X] -=1;
    }
  }         
} 



void YLOCATION(PARA_DATA *para, REAL **var, REAL *flag, REAL *y, REAL v0,
               int i, int j, int k,  REAL *OL, int *OC, int *LOC , int *COOD) {
  int imax = para->geom->imax, jmax = para->geom->jmax;
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  REAL *v=var[VY];
  REAL dt=para->mytime->dt;
  REAL *flagv=var[FLAGV];

  if(OL[Y]==y[FIX(OC[X],OC[Y],OC[Z])]) {      //velocity =0;
    COOD[Y]=0;
  }
  else if(OL[Y]<y[FIX(OC[X],OC[Y],OC[Z])])    //velocity >0
  {
    if(OC[Y]>0) OC[Y] -=1;                    //search backward
    if(OL[Y]>=y[FIX(OC[X],OC[Y],OC[Z])]) COOD[Y]=0; // comparing the coodinate
    if(flag[FIX(OC[X],OC[Y],OC[Z])]>0 || flagv[FIX(OC[X],OC[Y],OC[Z])]==4){ // hit boundary
      if(flagv[FIX(OC[X],OC[Y],OC[Z])]==0){  //inlet
        if(OL[Y]<y[FIX(OC[X],OC[Y],OC[Z])]) OL[Y]=y[FIX(OC[X],OC[Y],OC[Z])];
      }
      else {
        OL[Y]=y[FIX(OC[X],OC[Y]+1,OC[Z])]; // solid wall
      }
      LOC[Y]=0;
      COOD[Y]=0;
      OC[Y] +=1; 
    }
  }                             
  else {                        // velocity < 0
    if(OC[Y]<=jmax) OC[Y] +=1; //search in opposite direction
    if(OL[Y] <=y[FIX(OC[X],OC[Y],OC[Z])]) COOD[Y]=0; // comparing the coodinate
    if(flag[FIX(OC[X],OC[Y],OC[Z])]>0 || flagv[FIX(OC[X],OC[Y]-1,OC[Z])]==4){ // hit boundary
      if(flagv[FIX(OC[X],OC[Y]-1,OC[Z])]==0) { //inlet
        if(OL[Y]>y[FIX(OC[X],OC[Y],OC[Z])]) OL[Y]=y[FIX(OC[X],OC[Y],OC[Z])];
      }
      else {
        OL[Y]=y[FIX(OC[X],OC[Y]-1,OC[Z])];
        //LOC[Y]=0;
      }
      LOC[Y]=0;
      COOD[Y]=0;
      OC[Y] -=1; 
    }
  }
}


void ZLOCATION(PARA_DATA *para, REAL **var, REAL *flag, REAL *z, REAL w0,
               int i, int j, int k,  REAL *OL, int *OC, int *LOC , int *COOD){
  int imax = para->geom->imax, jmax = para->geom->jmax;
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  REAL *w=var[VZ];
  REAL dt=para->mytime->dt;
  REAL *flagw=var[FLAGW];

  if(OL[Z]==z[FIX(OC[X],OC[Y],OC[Z])]){
    COOD[Z]=0;
  }
  else if(OL[Z]<z[FIX(OC[X],OC[Y],OC[Z])])
  {
    if(OC[Z]>0) OC[Z] -=1;
    if(OL[Z]>=z[FIX(OC[X],OC[Y],OC[Z])]) COOD[Z]=0;   
    if(flag[FIX(OC[X],OC[Y],OC[Z])]>0 || flagw[FIX(OC[X],OC[Y],OC[Z])]==4 ) {
      if(flagw[FIX(OC[X],OC[Y],OC[Z])]==0){
        if(OL[Z]<z[FIX(OC[X],OC[Y],OC[Z])])  OL[Z]=z[FIX(OC[X],OC[Y],OC[Z])];
      }
      else {
        OL[Z]=z[FIX(OC[X],OC[Y],OC[Z]+1)];
      }
      LOC[Z]=0;
      COOD[Z]=0;
      OC[Z] +=1;
    }
  }
  else {
    if(OC[Z]<=kmax) OC[Z] +=1;
    if(OL[Z] <=z[FIX(OC[X],OC[Y],OC[Z])]) COOD[Z]=0;
    if(flag[FIX(OC[X],OC[Y],OC[Z])]>0 || flagw[FIX(OC[X],OC[Y],OC[Z]-1)]==4 ) { 
      if(flagw[FIX(OC[X],OC[Y],OC[Z]-1)]==0) {
        if(OL[Z]>z[FIX(OC[X],OC[Y],OC[Z])])  OL[Z]=z[FIX(OC[X],OC[Y],OC[Z])];
      }
      else {
        OL[Z]=z[FIX(OC[X],OC[Y],OC[Z]-1)]; 
      }
      LOC[Z]=0;
      COOD[Z]=0;
      OC[Z] -=1;
    }
  }
}



void XLOCATION_U(PARA_DATA *para, REAL **var, REAL *flag, REAL *x, REAL u0,
               int i, int j, int k,  REAL *OL, int *OC, int *LOC , int *COOD){

  int imax = para->geom->imax, jmax = para->geom->jmax; 
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  REAL *u=var[VX];
  REAL dt=para->mytime->dt;
  REAL *flagu=var[FLAGU];

  if(OL[X]==x[FIX(OC[X],OC[Y],OC[Z])]){
    COOD[X]=0;
  }
  else if(OL[X]<x[FIX(OC[X],OC[Y],OC[Z])])
  {
    if(OC[X]>0) OC[X] -=1;
    if(OL[X]>=x[FIX(OC[X],OC[Y],OC[Z])]) COOD[X]=0;
    if(flag[FIX(OC[X],OC[Y],OC[Z])]>=0){
      if(flag[FIX(OC[X],OC[Y],OC[Z])]==0){
        if(OL[X]<x[FIX(OC[X],OC[Y],OC[Z])]) OL[X]=x[FIX(OC[X],OC[Y],OC[Z])];
       }
      else {
        OL[X]=x[FIX(OC[X]+1,OC[Y],OC[Z])];
      }
      LOC[X]=0;
      COOD[X]=0;
      OC[X] +=1;
    }
  }
  else {
    if(OC[X]<=imax) OC[X] +=1;
    if(OL[X] <=x[FIX(OC[X],OC[Y],OC[Z])]) COOD[X]=0;
    if(flag[FIX(OC[X],OC[Y],OC[Z])]>=0){
      if(flagu[FIX(OC[X],OC[Y],OC[Z])]==0){
        if(OL[X]>x[FIX(OC[X],OC[Y],OC[Z])]) OL[X]=x[FIX(OC[X],OC[Y],OC[Z])];
      }
      else {
        OL[X]=x[FIX(OC[X]-1,OC[Y],OC[Z])];
      }
      LOC[X]=0;
      COOD[X]=0;
      OC[X] -=1;
    }
  }         
}



void YLOCATION_V(PARA_DATA *para, REAL **var, REAL *flag, REAL *y, REAL v0,
               int i, int j, int k,  REAL *OL, int *OC, int *LOC , int *COOD) {
  int imax = para->geom->imax, jmax = para->geom->jmax;
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  REAL *v=var[VY];
  REAL dt=para->mytime->dt;
  REAL *flagv=var[FLAGV];

  if(OL[Y]==y[FIX(OC[X],OC[Y],OC[Z])]) {      //velocity =0;
    COOD[Y]=0;
  }
  else if(OL[Y]<y[FIX(OC[X],OC[Y],OC[Z])])    //velocity >0
  {
    if(OC[Y]>0) OC[Y] -=1;                    //search backward
    if(OL[Y]>=y[FIX(OC[X],OC[Y],OC[Z])]) COOD[Y]=0; // comparing the coodinate
    if(flag[FIX(OC[X],OC[Y],OC[Z])]>=0){ // hit boundary
      if(flag[FIX(OC[X],OC[Y],OC[Z])]==0){  //inlet
        if(OL[Y]<y[FIX(OC[X],OC[Y],OC[Z])]) OL[Y]=y[FIX(OC[X],OC[Y],OC[Z])];
      }
      else {
        OL[Y]=y[FIX(OC[X],OC[Y]+1,OC[Z])]; // solid wall
      }
      LOC[Y]=0;
      COOD[Y]=0;
      OC[Y] +=1; 
    }
  }                             
  else {                        // velocity < 0
    if(OC[Y]<=jmax) OC[Y] +=1; //search in opposite direction
    if(OL[Y] <=y[FIX(OC[X],OC[Y],OC[Z])]) COOD[Y]=0; // comparing the coodinate
    if(flag[FIX(OC[X],OC[Y],OC[Z])]>=0){ // hit boundary
      if(flag[FIX(OC[X],OC[Y],OC[Z])]==0) { //inlet
        if(OL[Y]>y[FIX(OC[X],OC[Y],OC[Z])]) OL[Y]=y[FIX(OC[X],OC[Y],OC[Z])];
      }
      else {
        OL[Y]=y[FIX(OC[X],OC[Y]-1,OC[Z])];
      }
      LOC[Y]=0;
      COOD[Y]=0;
      OC[Y] -=1; 
    }
  }
}


void ZLOCATION_W(PARA_DATA *para, REAL **var, REAL *flag, REAL *z, REAL w0,
               int i, int j, int k,  REAL *OL, int *OC, int *LOC , int *COOD){
  int imax = para->geom->imax, jmax = para->geom->jmax;
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  REAL *w=var[VZ];
  REAL dt=para->mytime->dt;
  REAL *flagw=var[FLAGW];

  if(OL[Z]==z[FIX(OC[X],OC[Y],OC[Z])]){
    COOD[Z]=0;
  }
  else if(OL[Z]<z[FIX(OC[X],OC[Y],OC[Z])])
  {
    if(OC[Z]>0) OC[Z] -=1;
    if(OL[Z]>=z[FIX(OC[X],OC[Y],OC[Z])]) COOD[Z]=0;   
    if(flag[FIX(OC[X],OC[Y],OC[Z])]>=0 ) {
      if(flag[FIX(OC[X],OC[Y],OC[Z])]==0){
        if(OL[Z]<z[FIX(OC[X],OC[Y],OC[Z])])  OL[Z]=z[FIX(OC[X],OC[Y],OC[Z])];
      }
      else {
        OL[Z]=z[FIX(OC[X],OC[Y],OC[Z]+1)];
      }
      LOC[Z]=0;
      COOD[Z]=0;
      OC[Z] +=1;
    }
  }
  else {
    if(OC[Z]<=kmax) OC[Z] +=1;
    if(OL[Z] <=z[FIX(OC[X],OC[Y],OC[Z])]) COOD[Z]=0;
    if(flag[FIX(OC[X],OC[Y],OC[Z])]>=0 ) { 
      if(flagw[FIX(OC[X],OC[Y],OC[Z])]==0) {
        if(OL[Z]>z[FIX(OC[X],OC[Y],OC[Z])])  OL[Z]=z[FIX(OC[X],OC[Y],OC[Z])];
      }
      else {
        OL[Z]=z[FIX(OC[X],OC[Y],OC[Z]-1)]; 
      }
      LOC[Z]=0;
      COOD[Z]=0;
      OC[Z] -=1;
    }
  }
}
 
