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

FILE *file2;
FILE *file3;
FILE *file4;

/******************************************************************************
| FFD Solver
******************************************************************************/
void FFD_solver(PARA_DATA *para, REAL **var,int **BINDEX)
{
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
  REAL *T_mon,*U_mon,*P_mon;
  REAL T_mean,U_mean,P_mean;
  int i1=22,i2=22,j1=8,j2=8,k1=8,k2=8; // i1=3,i2=3,j1=1,j2=4,k1=5,k2=5; 

  int cal_mean = para->outp->cal_mean;


    file2 = fopen("velocity.txt", "w");
    if(file2 == NULL)
    {
    printf("Error:can not open input file!\n");  
    }
    else
     printf("file open successfully\n");

    file3 = fopen("temp.txt", "w");
    if(file3 == NULL)
    {
    printf("Error:can not open input file!\n");  
    }
    else
     printf("file open successfully\n");

    file4 = fopen("ventilation rate.txt", "w");
    if(file3 == NULL)
    {
    printf("Error:can not open input file!\n");  
    }
    else
     printf("file open successfully\n");
  
   T_mon =   (REAL *) malloc((2)*sizeof(REAL));
   U_mon =   (REAL *) malloc((2)*sizeof(REAL));
   P_mon =   (REAL *) malloc((2)*sizeof(REAL));
     if ( !T_mon ||!U_mon ||!P_mon) 
  {
    printf ("cannot allocate data\n" );
  }

  /*---------------------------------------------------------------------------
  | zero-equation length scale
  ---------------------------------------------------------------------------*/
   lengthscale(para,var);
  
  /*---------------------------------------------------------------------------
  | Solver Loop
  ---------------------------------------------------------------------------*/
  while( para->mytime->t_step < t_output)
  {
      vis_turb(para, var);      
      vel_step(para, var, BINDEX);
      temp_step(para, var,BINDEX);
      timing(para);
	  
      
    //U_mean=para->prob->resu;
    //P_mean=check_energy(para, var, BINDEX);
    //T_mean=check_avg(para,var,temp);

    //if((para->mytime->t_step%10)==0)
    //  {
    //            tsize= para->mytime->t_step/10;

    //     T_mon=(REAL *)realloc(T_mon,(tsize+1)*sizeof(REAL));
    //      U_mon=(REAL *)realloc(U_mon,(tsize+1)*sizeof(REAL));
    //      P_mon=(REAL *)realloc(P_mon,(tsize+1)*sizeof(REAL));

    //    if(!T_mon||!U_mon ||!P_mon )
    //           {
    //              printf("Reallocation failed\n");
    //              exit(1); 
    //      }

    //             T_mon[tsize]=T_mean;
    //             U_mon[tsize]=U_mean;
    //             P_mon[tsize]=P_mean;
    //   }
    //   
     
 

  } // End of While loop  

  /*---------------------------------------------------------------------------
  | Post Process
  ---------------------------------------------------------------------------*/
   cfl(para, var);
  //write_plume_mass(para, var, BINDEX);
  write_unsteady(para, var, "unsteady");
  write_Tecplot(para, var, "result",BINDEX);
  write_temperature(para, var,BINDEX);
  // write_time( tsize, T_mon,U_mon,P_mon);
   
 //  printf("iteru=%d\t iterp=%d\n",para->prob->iteru,para->prob->iterp);
    fclose(file2);
    fclose(file3);
    fclose(file4);
  free(T_mon);
  free(U_mon);

  para->prob->output = 1;
} // End of FFD_solver( ) 


/******************************************************************************
| calculate the temperature
******************************************************************************/ 
void temp_step(PARA_DATA *para, REAL **var,int **BINDEX)
{
  REAL *T = var[TEMP], *T0 = var[TMP1],*T1 = var[TMP2];            
  REAL *flagp = var[FLAGP];
  int i;
  int j=para->geom->jplume+1;
  int k=para->geom->kplume;
  int imax=para->geom->imax, jmax=para->geom->jmax;
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);

 set_bnd_temp(para, var, TEMP, T,BINDEX); 

 advection(para, var, TEMP,flagp, T0, T, BINDEX);

 if(para->prob->plume_mod==1) plume_thermal(para,var,BINDEX);

 //for(i=1;i<=imax;i++) {
 //  fprintf(file3,"%f\t",T0[FIX(i,j,k)]);
 //}
 //fprintf(file3,"\n");
 
 psi_conservation(para, var,T0,T,BINDEX);
 //for(i=1;i<=imax;i++) {
 //  fprintf(file3,"%f\t",T0[FIX(i,j,k)]);
 //}
 //fprintf(file3,"\n");

 diffusion(para, var, TEMP, T, T0, BINDEX);
 //printf("check energy=%f\n",check_avg(para, var, T)-check_avg(para, var, T0));
 //for(i=1;i<=imax;i++) {
 //  fprintf(file3,"%f\t",T[FIX(i,j,k)]);
 //}
 //fprintf(file3,"\n");
 //fprintf(file3,"\n");
 


 //printf("qheatsource=%f\n",qheatsource(para, var,BINDEX));
  //  plume(para,var,BINDEX);
 //source(para, var, TEMP, T, T0, BINDEX);

} // End of temp_step( )

/******************************************************************************
| calculate the temperature
******************************************************************************/ 
void den_step(PARA_DATA *para, REAL **var,int **BINDEX)
{
  REAL *den = var[DEN], *den0 = var[TMP1];
  REAL *flagp = var[FLAGP];
      
  set_bnd_density(para, var,DEN,den,BINDEX);
  advection(para, var, DEN,flagp, den0, den,BINDEX); 
  //psi_conservation(para, var,den0,den,BINDEX);

   //  swapT(para,var);

  diffusion(para, var, DEN, den, den0,BINDEX); 

   //   swapT(para,var);

} // End of den_step( )

/******************************************************************************
| calculate the velocity
******************************************************************************/ 

void vel_step(PARA_DATA *para, REAL **var,int **BINDEX)
{
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

  //for(i=1;i<=imax;i++) {
  //for(j=1;j<=jmax;j++) {
  //  fprintf(file2,"%f\t",v0[FIX(i,j,k)]);
  //}
  //fprintf(file2,"\n");

  diffusion(para, var, VX, u1, u0, BINDEX);   
  diffusion(para, var, VY, v1, v0, BINDEX); 
  diffusion(para, var, VZ, w1, w0, BINDEX); 

  //for(i=1;i<=imax;i++) {
  //for(j=1;j<=jmax;j++) {
  // fprintf(file2,"%f\t",v1[FIX(i,j,k)]);
  //}
  //fprintf(file2,"\n");

  if(para->bc->NBOUT!=0) mass_conservation(para, var,BINDEX);
  swapuvw(para,var);

  //if(para->prob->plume_mod==1) plume(para,var,BINDEX);
  //swapuvw_temp(para,var); 

  //for(j=1;j<=jmax;j++) {
  // fprintf(file2,"%f\t",v[FIX(i,j,k)]);
  //}
  //fprintf(file2,"\n");

  projection(para, var,BINDEX);

  if(para->prob->plume_mod==1) plume(para,var,BINDEX);

  //for(i=1;i<=imax;i++) {
  //for(j=1;j<=jmax;j++) {
  // fprintf(file2,"%f\t",v[FIX(i,j,k)]);
  //}
  //fprintf(file2,"\n");
  //fprintf(file2,"\n");

} // End of vel_step( )



/******************************************************************************
| Solver of equations
******************************************************************************/ 
void equ_solver(PARA_DATA *para, REAL **var, int var_type, REAL *psi) {
  REAL *flagp = var[FLAGP],*flagu = var[FLAGU],*flagv = var[FLAGV],*flagw = var[FLAGW];
  REAL res;

  switch (var_type) {
    case VX:
      Gauss_Seidel(para, var, flagu, psi);
      //res= check_residual(para,var,flagu,psi);
    break;
    case VY:
      Gauss_Seidel(para, var, flagv, psi);
      //res= check_residual(para,var,flagv,psi );
    break;
    case VZ:
      Gauss_Seidel(para, var, flagw, psi);
      // res= check_residual(para,var,flagw,psi );
    break;
    case TEMP:
    case IP:
    case DEN:
      GS_P(para, var, TEMP, psi);
      //res= check_residual(para,var,flagp,psi );
    break;
  }
     


}// end of equ_solver