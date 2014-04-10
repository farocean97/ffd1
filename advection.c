///////////////////////////////////////////////////////////////////////////////
///
/// \file   advection.c
///
/// \brief  Solver for advection step
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
/// This file provides functions that used for the advection step of FFD method.
/// The advection starts with \c advect(). Then different subroutines are 
/// called according to the properties of the variables that are sorted by
/// the location of variables assigned in the control volume. 
/// Velocities at X, Y and Z directions are locatted
/// on the surface of the control volume. They are computed using 
/// subroutines: \c traceback_uvw()
/// Scalar variables are in the center of control volume and they are computed
/// using \c traceback().
///
///////////////////////////////////////////////////////////////////////////////


#include <math.h>
#include <stdio.h>
#include "data_structure.h"
#include "boundary.h"
#include "advection.h"
#include "solver.h"
#include "utility.h"
#include "interpolation.h"

///////////////////////////////////////////////////////////////////////////////
///\brief  Entrance of advection step
///
/// Specific method for advection will be selected according to the variable 
/// type.
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param var_type The type of variable for advection solver
///\param flag Pointer to the cell flag
///\param d Pointer to the computed variables at previous time step
///\param d0 Pointer to the computed variables for current time step
///\param BINDEX Pointer to boundary index
///
///\return void No return needed
///////////////////////////////////////////////////////////////////////////////

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
      // For temperautre/Density,traceback from cell center
      traceback(para,var,BINDEX);

      FOR_EACH_CELL

        // flag=1 this is the solid cell
        if(flag[FIX(i,j,k)]==1) continue;
        // flag=0 this is the heat source cell
        else if(flag[FIX(i,j,k)]==0){
          // For heatsource cell,if plume model==1, plume model is applied
          if(para->prob->plume_mod==1) {
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
        // flag!=0& flag !=1 this is the fluid cell
        else
        {
          p=BINDEX[4][FIX(i,j,k)];
          q=BINDEX[5][FIX(i,j,k)];
          r=BINDEX[6][FIX(i,j,k)];

          d[FIX(i,j,k)] = interpolation_temp(para,var, d0, fx[FIX(i,j,k)], 
                                  fy[FIX(i,j,k)], fz[FIX(i,j,k)],i,j,k, p,q,r);
 
        }
      END_FOR
      break;

    // For U cell,traceback from cell surface
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
    // For V cell,traceback from cell surface
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
    // For W cell,traceback from cell surface
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
      
}// End of advection( )



///////////////////////////////////////////////////////////////////////////////
/// \brief Traceback for the variable at cell center
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param BINDEX Pointer to boundary index
///
///\return void No return needed
///////////////////////////////////////////////////////////////////////////////
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


  FOR_EACH_CELL
    if(flagp[FIX(i,j,k)]>0) continue;
    u0 = 0.5f*(u[FIX(i,j,k)]+u[FIX(i-1,j,k)]);
    v0 = 0.5f*(v[FIX(i,j,k)]+v[FIX(i,j-1,k)]);
    w0 = 0.5f*(w[FIX(i,j,k)]+w[FIX(i,j,k-1)]);
          
    //the location of the departure point
    OL[X] = x[FIX(i,j,k)] - u0*dt; 
    OL[Y] = y[FIX(i,j,k)] - v0*dt;
    OL[Z] = z[FIX(i,j,k)] - w0*dt;
    
    //Initial coordinate of the departure cell
    OC[X] = i; OC[Y] = j; OC[Z] = k;

    //COOD=1: the traceback in this direction is not completed.
    //COOD=0: the traceback in this direction is completed.
    COOD[X] =1; COOD[Y]=1; COOD[Z]=1;

    //LOC=1: the traceback doesn't cross boundary
    //LOC=0: the traceback crosses boundary
    LOC[X]  =1; LOC[Y] =1; LOC[Z] =1;

    while(COOD[X]==1 || COOD[Y] ==1 || COOD[Z] == 1){
      //traceback in X direction
      if(COOD[X]==1 && LOC[X]==1)
        XLOCATION(para, var, flagp, x, u0, i, j, k, OL,OC, LOC ,COOD);
      //traceback in Y direction
      if(COOD[Y]==1 && LOC[Y]==1)
        YLOCATION(para, var, flagp, y, v0, i, j, k, OL,OC, LOC ,COOD); 
      //traceback in Z direction
      if(COOD[Z]==1 && LOC[Z]==1)
        ZLOCATION(para, var, flagp, z, w0, i, j, k, OL,OC, LOC ,COOD); 
       
    } 

    //When velocity is positive and traceback crossed boundary, the identified
    //coordinate of the departure cell should be ajusted by -1.
    if(u0>=0 && LOC[X] == 0) OC[X] -=1;
	  if(v0>=0 && LOC[Y] == 0) OC[Y] -=1;
	  if(w0>=0 && LOC[Z] == 0) OC[Z] -=1;

    //When velocity is negtive and traceback didn't crossed boundary, the identified
    //coordinate of the departure cell should be ajusted by -1.
    if(u0<0 && LOC[X]>0) OC[X] -=1;
    if(v0<0 && LOC[Y]>0) OC[Y] -=1;
    if(w0<0 && LOC[Z]>0) OC[Z] -=1;

     //Store the coordinate of departure cell to BINDEX
    BINDEX[4][FIX(i,j,k)]=OC[X];
    BINDEX[5][FIX(i,j,k)]=OC[Y];
    BINDEX[6][FIX(i,j,k)]=OC[Z];

    //Calculate the distance ratio between departure point and surrounding cell centers.
    fx[FIX(i,j,k)] = (OL[X]- x[FIX(OC[X],OC[Y],OC[Z])])
                    /(x[FIX(OC[X]+1,OC[Y],OC[Z])]- x[FIX(OC[X],OC[Y],OC[Z])]); 
    fy[FIX(i,j,k)] = (OL[Y]- y[FIX(OC[X],OC[Y],OC[Z])])
                    /(y[FIX(OC[X],OC[Y]+1,OC[Z])]- y[FIX(OC[X],OC[Y],OC[Z])]);
    fz[FIX(i,j,k)] = (OL[Z]- z[FIX(OC[X],OC[Y],OC[Z])])
                    /(z[FIX(OC[X],OC[Y],OC[Z]+1)]- z[FIX(OC[X],OC[Y],OC[Z])]);

  
  END_FOR

}//End of traceback()

///////////////////////////////////////////////////////////////////////////////
///\brief Traceback for the variable at cell surfaces
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param var_type The type of variable for advection solver
///\param flag Pointer to the cell flag
///\param BINDEX Pointer to boundary index
///
///\return void No return needed
///////////////////////////////////////////////////////////////////////////////
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
    //traceback for U cell
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

    //traceback for V cell
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

    //traceback for W cell
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
}// End of traceback_UVW()


///////////////////////////////////////////////////////////////////////////////
/// Find the location and coordinates of departure points at x direction
///
/// Conducting backward tracing for the particle's X-location and 
/// corresponding coordinates at the previous time step.
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param flag Pointer to the property of the cell
///\param x Pointer to the x coordinates of the cell
///\param u0 X-velocity at time (t-1) in location x(t) 
///\param i I-index for cell at time t at x(t) 
///\param j J-index for cell at time t at x(t)
///\param k K-index for cell at time t at x(t)
///\param OL Pointer to the locations of particle at time (t-1)
///\param OC Pointer to the coordinates of particle at time (t-1)
///\param LOC Pointer to flags recording if tracing back hits the boundary
///\param COOD Pointer to record the status of tracing back process
///
///\return void No return needed
///////////////////////////////////////////////////////////////////////////////  
void XLOCATION(PARA_DATA *para, REAL **var, REAL *flag, REAL *x, REAL u0,
               int i, int j, int k,  REAL *OL, int *OC, int *LOC , int *COOD){

  int imax = para->geom->imax, jmax = para->geom->jmax; 
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  REAL *u=var[VX];
  REAL dt=para->mytime->dt;
  REAL *flagu=var[FLAGU];

  //if U=0 stop traceback
  if(OL[X]==x[FIX(OC[X],OC[Y],OC[Z])]){
    COOD[X]=0;
  }
  //if U>0,start searching cells on the left
  else if(OL[X]<x[FIX(OC[X],OC[Y],OC[Z])])
  {
    //if the cell index is larger than zero, keep searching.
    if(OC[X]>0) OC[X] -=1;

    //if the cell coordiate is smaller than that of departure point, 
    //stop searching and set COOD[X]=0
    if(OL[X]>=x[FIX(OC[X],OC[Y],OC[Z])]) COOD[X]=0;

    //if the cell flag shows the cell is solid(flag>0) or pass 
    //through wall partiions(flagu==4), stop searching in x direciton
    //set COOD=0,and assign LOC=0 
    //reset the cell index to the one next to the bounary cell so that trace
    //back in other directions can still continue.
    if(flag[FIX(OC[X],OC[Y],OC[Z])]>0 || flagu[FIX(OC[X],OC[Y],OC[Z])]==4){

      LOC[X]=0;
      COOD[X]=0;
      OC[X] +=1;

      //if the cell flag shows inlet and the traceback crosses the inlet cell,
      //reset the departure coordinate to that of the inlet cell 
      if(flagu[FIX(OC[X],OC[Y],OC[Z])]==0){
        if(OL[X]<x[FIX(OC[X],OC[Y],OC[Z])]) OL[X]=x[FIX(OC[X],OC[Y],OC[Z])];
      }

      //if the cell is not inlet, reset the departure coordinate to the one next
      //to the inlet cell
      else {
        OL[X]=x[FIX(OC[X]+1,OC[Y],OC[Z])];
      }
    }
  }

  //if the U<0, then search to the right direction.
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
 
