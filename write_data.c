/* ----------------------------------------------------------------------------
  
  Filename:      write.c

  Written by:   Wangda Zuo

  Task:          output the data as format for tecplot  

---------------------------------------------------------------------------- */  
#define _CRT_SECURE_NO_DEPRECATE
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "data_structure.h"
#include "write_data.h"
#include "utility.h"

FILE *file1;

/******************************************************************************
| Write the data to a file for Tecplot 
******************************************************************************/
int write_Tecplot(PARA_DATA *para, REAL **var, char *name,int **BINDEX)
{
  int i, j, k;
  int iter;
  int SI,SJ,SK,EI,EJ,EK;
  int imax=para->geom->imax, jmax=para->geom->jmax;
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  int n = para->mytime->t_step;
  int zone_num=para->geom->zone_wall;
  REAL *x = var[X], *y = var[Y], *z =var[Z];
  REAL *gx = var[GX], *gy = var[GY], *gz =var[GZ];
  REAL *aw = var[AW], *ae = var[AE], *as = var[AS], *an = var[AN];
  REAL *af = var[AF], *ab = var[AB],*ap=var[AP],*b=var[B];
  REAL *u0 = var[VX], *v0 = var[VY], *w0 = var[VZ], *p0 = var[IP], *d0 = var[DEN],*T0 = var[TEMP];
  REAL *u = var[VXM], *v = var[VYM], *w = var[VZM], *p= var[TMP1],*d = var[TMP2],*T = var[TEMPM];
  REAL *flagp = var[FLAGP],*flagu = var[FLAGU],*flagv = var[FLAGV],*flagw = var[FLAGW];
  REAL *vt=var[VT];


  char filename[20];  
  
  strcpy(filename, name);
  strcat(filename, ".plt");

  /* open output file */
  if((file1 = fopen( filename, "w" ))==NULL)
  {
    fprintf(stderr,"Error:can not open input file!\n");
    return -1;
  }

  corners(para, var,  u0);
  corners(para, var,  v0);
  corners(para, var,  w0);
  corners(para, var,  p0);
  corners(para, var,  d0);
  corners(para, var,  T0);

  u_vortex(para,var,u,u0);
  v_vortex(para,var,v,v0);
  w_vortex(para,var,w,w0);
  p_vortex(para,var,p,p0);
  p_vortex(para,var,T,T0);
  p_vortex(para,var,d,d0);


  corners_vortex(para, var,  u);
  corners_vortex(para, var,  v);
  corners_vortex(para, var,  w);
  corners_vortex(para, var,  p);
  corners_vortex(para, var,  d);
  corners_vortex(para, var,  T);


  uvw_center(para,var);

  
  fprintf( file1, "TITLE = ");

  fprintf( file1, "\"dt=%fs, t=%fs, nu=%f, Lx=%d, Ly=%d, Lz%d, Nx=%d, Ny=%d, Nz=%d \"\n",
           para->mytime->dt, para->mytime->t, para->prob->nu, para->geom->Lx, para->geom->Ly, para->geom->Lz,
           imax+w, jmax+2, kmax+2);

  fprintf( file1, 
           "VARIABLES =X, Y, Z, U, V, W,P,T,DEN \n");
  fprintf( file1, "ZONE F=POINT, I=%d, J=%d, K=%d\n", imax+2, jmax+2, kmax+2 );

 
    for(k=0; k<=kmax+1; k++) 
      for(j=0; j<=jmax+1; j++) 
          for(i=0; i<=imax+1; i++) 
         {

               fprintf( file1, "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", x[FIX(i,j,k)], y[FIX(i,j,k)],
               z[FIX(i,j,k)],u0[FIX(i,j,k)],v0[FIX(i,j,k)],w0[FIX(i,j,k)],p0[FIX(i,j,k)],T0[FIX(i,j,k)],var[DENS][FIX(i,j,k)]);   
        }



   for(iter=1;iter<=zone_num;iter++)
      {
       SI=BINDEX[7][iter];
       EI=BINDEX[8][iter];
       SJ=BINDEX[9][iter];
       EJ=BINDEX[10][iter];
       SK=BINDEX[11][iter];
       EK=BINDEX[12][iter];

       fprintf (file1,"ZONE T= \"Subzone%d \"\n",iter);
             fprintf( file1, "I=%d, J=%d, K=%d,ZONETYPE=Ordered\n", EI-SI+1, EJ-SJ+1, EK-SK+1);
             fprintf( file1, "DATAPACKING=POINT\n");

       for(k=SK; k<=EK; k++) 
           for(j=SJ; j<=EJ; j++)
           for(i=SI; i<=EI; i++)
           {
                 fprintf( file1, "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", gx[FIX(i,j,k)], gy[FIX(i,j,k)],
                 gz[FIX(i,j,k)],u[FIX(i,j,k)],v[FIX(i,j,k)],w[FIX(i,j,k)],p[FIX(i,j,k)],T[FIX(i,j,k)],d[FIX(i,j,k)]);    
          }
       
       
       }
       
     
 
  fclose(file1);

  printf("The data file %s has been written!\n", name);
  return 0;

} //write_data()


int write_data1(PARA_DATA *para, REAL **var, char *name)
{
  int i, j, k;
  int imax=para->geom->imax, jmax=para->geom->jmax;
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  int n = para->mytime->t_step;
  REAL *x = var[X], *y = var[Y], *z =var[Z];
  REAL *u = var[VXM], *v = var[VYM], *w = var[VZM], *d = var[DEN];
  REAL *T = var[TEMP], *Tm = var[TEMPM];
  REAL *tmp1 = var[TMP1], *tmp2 = var[TMP2], *tmp3 = var[TMP3];
  REAL *flagp = var[FLAGP],*flagu = var[FLAGU],*flagv = var[FLAGV],*flagw = var[FLAGW];
  REAL *aw = var[AW], *ae = var[AE], *as = var[AS], *an = var[AN];
  REAL *af = var[AF], *ab = var[AB],*ap0=var[AP0], *b=var[B],*ap=var[AP];
  char filename[20];  
  
  strcpy(filename, name);
  strcat(filename, ".plt");

  /* open output file */
  if((file1 = fopen( filename, "w" ))==NULL)
  {
    fprintf(stderr,"Error:can not open input file!\n");
    return -1;
  }

  /*  for(j=0; j<=jmax+1; j++)
	  {
		for(k=0; k<=kmax+1; k++)
          {
            
                u[FIX(imax+1,j,k)] = u[FIX(imax,j,k)];

                     for(i=imax; i>=1; i--)
                     {
                          u[FIX(i,j,k)] = 0.5f * (u[FIX(i,j,k)]+u[FIX(i-1,j,k)]);
                      }
            }
	  }


	
    for(i=0; i<=imax+1; i++)
	  {
		for(k=0; k<=kmax+1; k++)
          {
            
                v[FIX(i,jmax+1,k)] = v[FIX(i,jmax,k)];

                     for(j=jmax; j>=1; j--)
                     {
                          v[FIX(i,j,k)] = 0.5f * (v[FIX(i,j,k)]+v[FIX(i,j-1,k)]);
                        }
            }
	  }
   

	 for(i=0; i<=imax+1; i++)
	  {
		for(j=0; j<=jmax+1; j++)
          {
            
                w[FIX(i,j,kmax+1)] = w[FIX(i,j,kmax)];

                     for(k=kmax; k>=1; k--)
                     {
                          w[FIX(i,j,k)] = 0.5f * (w[FIX(i,j,k)]+w[FIX(i,j,k-1)]);
                      }
            }
	  }
  */   
  fprintf( file1, "TITLE = ");

  fprintf( file1, "\"dt=%fs, t=%fs, nu=%f, Lx=%d, Ly=%d, Lz%d, Nx=%d, Ny=%d, Nz=%d \"\n",
           para->mytime->dt, para->mytime->t, para->prob->nu, para->geom->Lx, para->geom->Ly, para->geom->Lz,
           imax+2, jmax+2, kmax+2);

  fprintf( file1, 
           "VARIABLES =X, Y, Z, U, V, W,AN,AB,AF \n");
  fprintf( file1, "ZONE F=POINT, I=%d, J=%d, K=%d\n", imax+2, jmax+2, kmax+2 );

   FOR_ALL_CELL

    fprintf( file1, "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",
      x[FIX(i,j,k)],y[FIX(i,j,k)],z[FIX(i,j,k)],u[FIX(i,j,k)], v[FIX(i,j,k)], w[FIX(i,j,k)], an[FIX(i,j,k)], ab[FIX(i,j,k)], b[FIX(i,j,k)]);   
  END_FOR
  

  fclose(file1);
  
  printf("The data file %s has been written!\n", name);
  return 0;

} //write_data()

int write_unsteady(PARA_DATA *para, REAL **var, char *name)
{
  int i,j,k;
  int imax=para->geom->imax, jmax=para->geom->jmax;
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  REAL *u = var[VX], *v = var[VY], *w = var[VZ], *p = var[IP];
  REAL  *d = var[DEN];
  REAL *T = var[TEMP];

  char filename[20];  
  
  strcpy(filename, name);
  strcat(filename, ".plt");

  /* open output file */
  if((file1 = fopen( filename, "w" ))==NULL)
  {
    fprintf(stderr,"Error:can not open input file!\n");
    return -1;
  }
 
  FOR_ALL_CELL
    fprintf( file1, "%f\t%f\t%f\t%f\t%f\t%f\n",u[FIX(i,j,k)], v[FIX(i,j,k)], w[FIX(i,j,k)], T[FIX(i,j,k)],d[FIX(i,j,k)], p[FIX(i,j,k)]);    
  END_FOR
  
  fclose(file1);
  printf("The data file %s has been written!\n", name);
  return 0;

} //write_data()

int write_SCI(PARA_DATA *para, REAL **var, char *name)
{
  int i, j, k;
  int IPR,IU,IV,IW,IT,IC1,IC2,IC3,IC4,IC5,IC6,IC7;
  int imax=para->geom->imax, jmax=para->geom->jmax;
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  int n = para->mytime->t_step;
  REAL *x = var[X], *y = var[Y], *z =var[Z];
  REAL *u = var[VX], *v = var[VY], *w = var[VZ], *p = var[IP];
  REAL *um = var[VXM], *vm = var[VYM], *wm = var[VZM], *d = var[DEN];
  REAL *T = var[TEMP], *Tm = var[TEMPM];
  REAL *tmp1 = var[TMP1], *tmp2 = var[TMP2], *tmp3 = var[TMP3];
  char filename[20];  

  
  strcpy(filename, name);
  strcat(filename, ".cfd");

  if((file1 = fopen( filename, "w" ))==NULL)
  {
    fprintf(stderr,"Error:can not open input file!\n");
    return -1;
  }

  IPR=1;
  IU=1;
  IV=1;
  IW=1;
  IT=1;
  IC1=0;
  IC2=0;
  IC3=0;
  IC4=0;
  IC5=0;
  IC6=0;
  IC7=0;


    for(j=0; j<=jmax+1; j++)
    {
    for(k=0; k<=kmax+1; k++)
          {
            
                u[FIX(imax+1,j,k)] = u[FIX(imax,j,k)];
                um[FIX(imax+1,j,k)] = um[FIX(imax,j,k)];   

                     for(i=imax; i>=1; i--)
                     {
                          u[FIX(i,j,k)] = 0.5f * (u[FIX(i,j,k)]+u[FIX(i-1,j,k)]);
                          um[FIX(i,j,k)] = 0.5f * (um[FIX(i,j,k)]+um[FIX(i-1,j,k)]);
                      }
            }
    }


  
    for(i=0; i<=imax+1; i++)
    {
    for(k=0; k<=kmax+1; k++)
          {
            
                v[FIX(i,jmax+1,k)] = v[FIX(i,jmax,k)];
                vm[FIX(i,jmax+1,k)] = vm[FIX(i,jmax,k)];  

                     for(j=jmax; j>=1; j--)
                     {
                          v[FIX(i,j,k)] = 0.5f * (v[FIX(i,j,k)]+v[FIX(i,j-1,k)]);
                          vm[FIX(i,j,k)] = 0.5f * (vm[FIX(i,j,k)]+vm[FIX(i,j-1,k)]);
                      }
            }
    }
   

   for(i=0; i<=imax+1; i++)
    {
    for(j=0; j<=jmax+1; j++)
          {
            
                w[FIX(i,j,kmax+1)] = w[FIX(i,j,kmax)];
                wm[FIX(i,j,kmax+1)] = wm[FIX(i,j,kmax)];  

                     for(k=kmax; k>=1; k--)
                     {
                          w[FIX(i,j,k)] = 0.5f * (w[FIX(i,j,k)]+w[FIX(i,j,k-1)]);
                          wm[FIX(i,j,k)] = 0.5f * (wm[FIX(i,j,k)]+wm[FIX(i,j,k-1)]);
                      }
            }
    }
  
 //W-S-B
  p[FIX(0,0,0)] = (p[FIX(0,1,0)]+p[FIX(1,0,0)]+p[FIX(0,0,1)]) / 3.0f;
  //W-N-B
  p[FIX(0,jmax+1,0)] = ( p[FIX(1,jmax+1,0)]+p[FIX(0,jmax,0)]
                         +p[FIX(0,jmax+1,1)]) / 3.0f;
  //E-S-B
  p[FIX(imax+1,0,0)] = ( p[FIX(imax,0,0)]+p[FIX(imax+1,1,0)]
                         +p[FIX(imax+1,0,1)]) / 3.0f;
  //E-N-B
  p[FIX(imax+1,jmax+1,0)] = ( p[FIX(imax,jmax+1,0)]+p[FIX(imax+1,jmax,0)]
                              +p[FIX(imax+1,jmax+1,1)]) / 3.0f;
  //W-S-F
  p[FIX(0,0,kmax+1)] = ( p[FIX(0,1,kmax+1)]+p[FIX(1,0,kmax+1)]
                         +p[FIX(0,0,kmax)]) / 3.0f;  
  //W-N-F
  p[FIX(0,jmax+1,kmax+1)] = ( p[FIX(1,jmax+1,kmax+1)]+p[FIX(0,jmax,kmax+1)]
                              +p[FIX(0,jmax+1,kmax)]) / 3.0f;

  //E-S-F
  p[FIX(imax+1,0,kmax+1)] = ( p[FIX(imax,0,kmax+1)]+p[FIX(imax+1,1,kmax+1)]
                              +p[FIX(imax+1,0,kmax)]) / 3.0f;
  //E-N-F
  p[FIX(imax+1,jmax+1,kmax+1)] = ( p[FIX(imax,jmax+1,0)]+p[FIX(imax+1,jmax,0)]
                                   +p[FIX(imax+1,jmax+1,kmax)]) / 3.0f;




  fprintf( file1, "%e\t%e\t%e\n", para->geom->Lx, para->geom->Ly, para->geom->Lz);
  fprintf( file1, "%d\t%d\t%d\n", imax,jmax,kmax);
  fprintf( file1, "%d\t%d\t%d\t%d\t%d\t%d\n", IPR,IU,IV,IW,IT,IC1);
  fprintf( file1, "%d\t%d\t%d\t%d\t%d\t%d\n", IC2,IC3,IC4,IC5,IC6,IC7);

  for(i=1;i<=imax;i++)  fprintf( file1, "%e\t", x[FIX(i,j,k)]);
  fprintf( file1, "\n");
  for(j=1;j<=jmax;j++)  fprintf( file1, "%e\t", y[FIX(i,j,k)]);
  fprintf( file1, "\n");
  for(k=1;k<=kmax;k++)  fprintf( file1, "%e\t", z[FIX(i,j,k)]);
  fprintf( file1, "\n");

  for(j=1;j<=jmax;j++)
    for(i=1;i<=imax;i++)
    {
        fprintf( file1, "%d\t%d\n", i,j);
   if(IPR==1)
   {
   for(k=1;k<=kmax;k++)  fprintf( file1, "%e\t", p[FIX(i,j,k)]);   
     fprintf( file1, "\n");
   }
   if(IU==1)
   {
   for(k=1;k<=kmax;k++)  fprintf( file1, "%e\t", u[FIX(i,j,k)]);   
     fprintf( file1, "\n");
   }
   if(IV==1)
   {
   for(k=1;k<=kmax;k++)  fprintf( file1, "%e\t", v[FIX(i,j,k)]);   
     fprintf( file1, "\n");
   }
   if(IW==1)
   {
   for(k=1;k<=kmax;k++)  fprintf( file1, "%e\t", w[FIX(i,j,k)]);   
     fprintf( file1, "\n");
   }

   if(IT==1)
   {
   for(k=1;k<=kmax;k++)  fprintf( file1, "%e\t", T[FIX(i,j,k)]);   
     fprintf( file1, "\n");
   }

   for(k=1;k<=kmax;k++)  fprintf( file1, "%e\t", T[FIX(i,j,k)]);  //turbulence intensity 
     fprintf( file1, "\n");
    }


  fclose(file1);
  
  printf("The data file %s has been written!\n", name);
  return 0;

} //write_data()

int write_time(int tsize, REAL *T_mon, REAL *U_mon,REAL *P_mon)
{
  int i;
  file1 = fopen("time.txt", "w");
    if(file1 == NULL)
    {
    printf("Error:can not open input file!\n");  /*fprintf(stderr,"Error:can not open input file!\n")*/
  return -1;
    }
    else
  printf("file open successfully\n");

 
   for (i=0;i<=tsize;i++)
   {
     fprintf (file1,"%d\t%f\t%.3g\t%.3g\n",i,T_mon[i],U_mon[i],P_mon[i]);
   }

   fclose(file1);
  
  printf("The data file result.plt has been written!\n");

  return 0;

} //write_data()

int write_displacement(PARA_DATA *para, REAL **var)
{
     int i, j, k;
     int imax=para->geom->imax, jmax=para->geom->jmax;
     int kmax = para->geom->kmax;
     int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
     REAL *x = var[X], *y = var[Y], *z =var[Z];
   REAL dy,dz;
   REAL integdist,integden,disp;
     REAL *den=var[DEN];

  file1 = fopen("displacement.txt", "w");
    if(file1 == NULL)
    {
    printf("Error:can not open input file!\n");  /*fprintf(stderr,"Error:can not open input file!\n")*/
  return -1;
    }
    else
  printf("file open successfully\n");

 
   for (i=0;i<=imax;i++)
   {
     integdist=0;
     integden=0;
       for(j=1;j<=jmax;j++)
       for(k=1;k<=kmax;k++)
       {
         dy=y[FIX(i,j,k)]-y[FIX(i,11,k)];
       dz=z[FIX(i,j,k)]-z[FIX(i,j,11)];
       integdist +=den[FIX(i,j,k)]*(dy*dy+dz*dz);
       integden +=den[FIX(i,j,k)];

       }

        disp=integdist/integden;
       
    fprintf (file1,"%d\t%.3g\t%.3g\n",i,x[FIX(i,j,k)],disp);
    }

   fclose(file1);
  
  printf("The data file result.plt has been written!\n");

  return 0;

} //write_data()



int write_data_block(PARA_DATA *para, REAL **var, char *name)
{
  int i, j, k;
  int index;
  int imax=para->geom->imax, jmax=para->geom->jmax;
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  int n = para->mytime->t_step;
  REAL *x = var[X], *y = var[Y], *z =var[Z];
  REAL *gx = var[GX], *gy = var[GY], *gz =var[GZ];
  REAL *aw = var[AW], *ae = var[AE], *as = var[AS], *an = var[AN];
  REAL *af = var[AF], *ab = var[AB],*ap=var[AP],*b=var[B];
  REAL *u = var[VX], *v = var[VY], *w = var[VZ], *p = var[IP];
  REAL *um = var[VXM], *vm = var[VYM], *wm = var[VZM], *d = var[DEN];
  REAL *T = var[TEMP], *Tm = var[TEMPM];
  REAL *flagp = var[FLAGP],*flagu = var[FLAGU],*flagv = var[FLAGV],*flagw = var[FLAGU];


  char filename[20];  
  
  strcpy(filename, name);
  strcat(filename, ".plt");

  /* open output file */
  if((file1 = fopen( filename, "w" ))==NULL)
  {
    fprintf(stderr,"Error:can not open input file!\n");
    return -1;
  }

    for(j=0; j<=jmax+1; j++)
    {
    for(k=0; k<=kmax+1; k++)
          {
            
                u[FIX(imax+1,j,k)] = u[FIX(imax,j,k)];

                     for(i=imax; i>=1; i--)
                     {
                          u[FIX(i,j,k)] = 0.5f * (u[FIX(i,j,k)]+u[FIX(i-1,j,k)]);
                      }
            }
    }


  
    for(i=0; i<=imax+1; i++)
    {
    for(k=0; k<=kmax+1; k++)
          {
            
                v[FIX(i,jmax+1,k)] = v[FIX(i,jmax,k)];

                     for(j=jmax; j>=1; j--)
                     {
                          v[FIX(i,j,k)] = 0.5f * (v[FIX(i,j,k)]+v[FIX(i,j-1,k)]);
                        }
            }
    }
   

   for(i=0; i<=imax+1; i++)
    {
    for(j=0; j<=jmax+1; j++)
          {
            
                w[FIX(i,j,kmax+1)] = w[FIX(i,j,kmax)];

                     for(k=kmax; k>=1; k--)
                     {
                          w[FIX(i,j,k)] = 0.5f * (w[FIX(i,j,k)]+w[FIX(i,j,k-1)]);
                      }
            }
    }
  
  //FOR_EACH_CELL
  //  if(i>i1 && i<=i2 && j>j1 && j<=j2 && k<=k2) T[FIX(i,j,k)]=308.3;
  //END_FOR


 //W-S-B
  p[FIX(0,0,0)] = (p[FIX(0,1,0)]+p[FIX(1,0,0)]+p[FIX(0,0,1)]) / 3.0f;
  //W-N-B
  p[FIX(0,jmax+1,0)] = ( p[FIX(1,jmax+1,0)]+p[FIX(0,jmax,0)]
                         +p[FIX(0,jmax+1,1)]) / 3.0f;
  //E-S-B
  p[FIX(imax+1,0,0)] = ( p[FIX(imax,0,0)]+p[FIX(imax+1,1,0)]
                         +p[FIX(imax+1,0,1)]) / 3.0f;
  //E-N-B
  p[FIX(imax+1,jmax+1,0)] = ( p[FIX(imax,jmax+1,0)]+p[FIX(imax+1,jmax,0)]
                              +p[FIX(imax+1,jmax+1,1)]) / 3.0f;
  //W-S-F
  p[FIX(0,0,kmax+1)] = ( p[FIX(0,1,kmax+1)]+p[FIX(1,0,kmax+1)]
                         +p[FIX(0,0,kmax)]) / 3.0f;  
  //W-N-F
  p[FIX(0,jmax+1,kmax+1)] = ( p[FIX(1,jmax+1,kmax+1)]+p[FIX(0,jmax,kmax+1)]
                              +p[FIX(0,jmax+1,kmax)]) / 3.0f;

  //E-S-F
  p[FIX(imax+1,0,kmax+1)] = ( p[FIX(imax,0,kmax+1)]+p[FIX(imax+1,1,kmax+1)]
                              +p[FIX(imax+1,0,kmax)]) / 3.0f;
  //E-N-F
  p[FIX(imax+1,jmax+1,kmax+1)] = ( p[FIX(imax,jmax+1,0)]+p[FIX(imax+1,jmax,0)]
                                   +p[FIX(imax+1,jmax+1,kmax)]) / 3.0f;

                    //W-S-B
  T[FIX(0,0,0)] = (T[FIX(0,1,0)]+T[FIX(1,0,0)]+T[FIX(0,0,1)]) / 3.0f;
  //W-N-B
  T[FIX(0,jmax+1,0)] = ( T[FIX(1,jmax+1,0)]+T[FIX(0,jmax,0)]
                         +T[FIX(0,jmax+1,1)]) / 3.0f;
  //E-S-B
  T[FIX(imax+1,0,0)] = ( T[FIX(imax,0,0)]+T[FIX(imax+1,1,0)]
                         +T[FIX(imax+1,0,1)]) / 3.0f;
  //E-N-B
  T[FIX(imax+1,jmax+1,0)] = ( T[FIX(imax,jmax+1,0)]+T[FIX(imax+1,jmax,0)]
                              +T[FIX(imax+1,jmax+1,1)]) / 3.0f;
  //W-S-F
  T[FIX(0,0,kmax+1)] = ( T[FIX(0,1,kmax+1)]+T[FIX(1,0,kmax+1)]
                         +T[FIX(0,0,kmax)]) / 3.0f;  
  //W-N-F
  T[FIX(0,jmax+1,kmax+1)] = ( T[FIX(1,jmax+1,kmax+1)]+T[FIX(0,jmax,kmax+1)]
                              +T[FIX(0,jmax+1,kmax)]) / 3.0f;

  //E-S-F
  T[FIX(imax+1,0,kmax+1)] = ( T[FIX(imax,0,kmax+1)]+T[FIX(imax+1,1,kmax+1)]
                              +T[FIX(imax+1,0,kmax)]) / 3.0f;
  //E-N-F
  T[FIX(imax+1,jmax+1,kmax+1)] = ( T[FIX(imax,jmax+1,0)]+T[FIX(imax+1,jmax,0)]
                                   +T[FIX(imax+1,jmax+1,kmax)]) / 3.0f;

  
  fprintf( file1, "TITLE = ");

  fprintf( file1, "\"dt=%fs, t=%fs, nu=%f, Lx=%d, Ly=%d, Lz%d, Nx=%d, Ny=%d, Nz=%d \"\n",
           para->mytime->dt, para->mytime->t, para->prob->nu, para->geom->Lx, para->geom->Ly, para->geom->Lz,
           imax+2, jmax+2, kmax+2);
  
  fprintf( file1, 
           "VARIABLES = X, Y, Z, U, V, W\n"); //

  fprintf (file1,"ZONE T= \"ZONE 001 \"\n");
  fprintf( file1, "I=%d, J=%d, K=%d,ZONETYPE=Ordered\n", imax+2, jmax+2, kmax+2 );
  fprintf( file1, "DATAPACKING=BLOCK\n");
  fprintf( file1, "VARLOCATION=([4-6]=CELLCENTERED)\n");

      FOR_ALL_CELL
         index=FIX(i,j,k)+1;
       fprintf( file1, "%f\t", var[GX][FIX(i,j,k)]);
       if(index%10==0) fprintf( file1, "\n");
    END_FOR
       if(index%10!=0) fprintf( file1, "\n");

      FOR_ALL_CELL
         index=FIX(i,j,k)+1;
       fprintf( file1, "%f\t", var[GY][FIX(i,j,k)]);
       if(index%10==0) fprintf( file1, "\n");
    END_FOR
       if(index%10!=0) fprintf( file1, "\n");

      FOR_ALL_CELL
         index=FIX(i,j,k)+1;
       fprintf( file1, "%f\t", var[GZ][FIX(i,j,k)]);
       if(index%10==0) fprintf( file1, "\n");
    END_FOR
       if(index%10!=0) fprintf( file1, "\n");

      FOR_ALL_CELL
         index=FIX(i,j,k)+1;
       fprintf( file1, "%f\t", var[VX][FIX(i,j,k)]);
       if(index%10==0) fprintf( file1, "\n");
    END_FOR
       if(index%10!=0) fprintf( file1, "\n");

      FOR_ALL_CELL
         index=FIX(i,j,k)+1;
       fprintf( file1, "%f\t", var[VY][FIX(i,j,k)]);
       if(index%10==0) fprintf( file1, "\n");
    END_FOR
       if(index%10!=0) fprintf( file1, "\n");

    FOR_ALL_CELL
         index=FIX(i,j,k)+1;
       fprintf( file1, "%f\t", var[VX][FIX(i,j,k)]);
       if(index%10==0) fprintf( file1, "\n");
    END_FOR
       if(index%10!=0) fprintf( file1, "\n");

  

  fclose(file1);

  printf("The data file %s has been written!\n", name);
  return 0;

} //write_data()

void write_temperature(PARA_DATA *para, REAL **var, int **BINDEX)
{
  int i, j, k;
  int imax=para->geom->imax, jmax=para->geom->jmax;
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  REAL *T=var[TEMP],*gx=var[GX],*gy=var[GY],*gz=var[GZ];
  REAL *y=var[Y];
  REAL *v=var[VY];
  REAL *flagp=var[FLAGP];
  REAL tmp,tmp_vis,tmp_loc,T_avg,TLOC_avg;
  REAL area,sum_area;
  REAL Vis_avg;
  REAL mass=0;
  int i_center=11;
  int k_center=11;

  file1 = fopen("average_temp.txt", "w");
    if(file1 == NULL)
    {
    printf("Error:can not open input file!\n");  /*fprintf(stderr,"Error:can not open input file!\n")*/
    }
    else
  printf("file open successfully\n");

    for(j=1;j<=jmax;j++) {
     tmp=0;tmp_loc=0;
      tmp_vis=0;
      sum_area=0;
      for(i=1;i<=imax;i++)
        for(k=1;k<=kmax;k++) {
          if(flagp[FIX(i,j,k)]==1 || BINDEX[20][FIX(i,j,k)]==1) continue;
          area= (gx[FIX(i,j,k)]-gx[FIX(i-1,j,k)])*(gz[FIX(i,j,k)]-gz[FIX(i,j,k-1)]);
          tmp += area*T[FIX(i,j,k)];
          sum_area += area;
        }
      T_avg = tmp/sum_area;

      
      sum_area=0;
     
     Vis_avg=v[FIX(i_center,j,k_center)];
     TLOC_avg = T[FIX(i_center,j,k_center)];

     mass=0;

     fprintf (file1,"%f\t%f\t%f\t%f\t%f\n",y[FIX(i_center,j,k_center)],Vis_avg,TLOC_avg,T_avg,mass);
    }

    fclose(file1);
  
    printf("The data file result.plt has been written!\n");


} //write_data()

void write_plume_mass(PARA_DATA *para, REAL **var, int **BINDEX) {
  int i, j, k;
  int imax = para->geom->imax, jmax = para->geom->jmax;
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  REAL *x = var[X], *y = var[Y], *z = var[Z];
  REAL *gx = var[GX], *gy = var[GY], *gz = var[GZ]; 
  REAL *u = var[VX], *v = var[VY], *w = var[VZ];  
  REAL *u1 = var[VXC], *v1 = var[VYC], *w1 = var[VZC];  
  REAL Dx,Dy,Dz;
  REAL b1,b2;
  int jplume;
  REAL zv=para->geom->zv;

    file1 = fopen("plume_mass.txt", "w");
    if(file1 == NULL)
    {
    printf("Error:can not open input file!\n");  
    }
    else
     printf("file open successfully\n");
  
    jplume=para->geom->jplume;
    i=para->geom->iplume;
    k=para->geom->kplume;

    for(j=jplume;j<=jmax;j++) {
    Dx = gx[FIX(i,j,k)]-gx[FIX(i-1,j,k)];
    Dy = gy[FIX(i,j ,k)]-gy[FIX(i,j-1,k)];
    Dz = gz[FIX(i,j ,k)]-gz[FIX(i,j,k-1)];
 
    if(BINDEX[20][FIX(i,j,k)]==1) {
        b1=6/5*0.09*(gy[FIX(i,j-1,k)]-zv);
        b2=6/5*0.09*(gy[FIX(i,j,k)]-zv);
        fprintf(file1,"%f\t",u[FIX(i-1,j,k)]);
        fprintf(file1,"%f\t",u[FIX(i,j,k)]);
        fprintf(file1,"%f\t",w[FIX(i,j,k-1)]);
        fprintf(file1,"%f\t",w[FIX(i,j,k)]);
        fprintf(file1,"%f\t",v[FIX(i,j-1,k)]);
        fprintf(file1,"%f\t",v[FIX(i,j,k)]);
        //fprintf(file1,"%f\t",3.14f*b1*b1*(1-exp(-Dx*Dz/(b1*b1)))*v1[FIX(i,j-1,k)]);
        //fprintf(file1,"%f\t",3.14f*b2*b2*(1-exp(-Dx*Dz/(b2*b2)))*v1[FIX(i,j,k)]);
        fprintf(file1,"\n");

        fprintf(file1,"%f\t",Dy*Dz*u[FIX(i-1,j,k)]);
        fprintf(file1,"%f\t",Dy*Dz*u[FIX(i,j,k)]);
        fprintf(file1,"%f\t",Dx*Dy*w[FIX(i,j,k-1)]);
        fprintf(file1,"%f\t",Dx*Dy*w[FIX(i,j,k)]);
        fprintf(file1,"%f\t",Dx*Dz*v[FIX(i,j-1,k)]);
        fprintf(file1,"%f\t",Dx*Dz*v[FIX(i,j,k)]);
        //fprintf(file1,"%f\t",3.14f*b1*b1*(1-exp(-Dx*Dz/(b1*b1)))*v[FIX(i,j-1,k)]);
        //fprintf(file1,"%f\t",3.14f*b2*b2*(1-exp(-Dx*Dz/(b2*b2)))*v[FIX(i,j,k)]);
         fprintf(file1,"\n");
         fprintf(file1,"\n");
      }
    else if(BINDEX[20][FIX(i,j-1,k)]==1) {
        b1=6/5*0.09*(gy[FIX(i,j-1,k)]-zv);
        b2=6/5*0.09*(gy[FIX(i,j,k)]-zv);
        fprintf(file1,"%f\t",u[FIX(i-1,j,k)]);
        fprintf(file1,"%f\t",u[FIX(i,j,k)]);
        fprintf(file1,"%f\t",w[FIX(i,j,k-1)]);
        fprintf(file1,"%f\t",w[FIX(i,j,k)]);
        fprintf(file1,"%f\t",v[FIX(i,j-1,k)]);
        fprintf(file1,"%f\t",v[FIX(i,j,k)]);
        fprintf(file1,"\n");

        fprintf(file1,"%f\t",Dy*Dz*u[FIX(i-1,j,k)]);
        fprintf(file1,"%f\t",Dy*Dz*u[FIX(i,j,k)]);
        fprintf(file1,"%f\t",Dx*Dy*w[FIX(i,j,k-1)]);
        fprintf(file1,"%f\t",Dx*Dy*w[FIX(i,j,k)]);
       // fprintf(file1,"%f\t",3.14f*b1*b1*(1-exp(-Dx*Dz/(b1*b1)))*v[FIX(i,j-1,k)]);
        fprintf(file1,"%f\t",Dx*Dz*v[FIX(i,j-1,k)]);
        fprintf(file1,"%f\t",Dx*Dz*v[FIX(i,j,k)]);
        fprintf(file1,"\n");
        fprintf(file1,"\n");
    }

    else {
        fprintf(file1,"%f\t",u[FIX(i-1,j,k)]);
        fprintf(file1,"%f\t",u[FIX(i,j,k)]);
        fprintf(file1,"%f\t",w[FIX(i,j,k-1)]);
        fprintf(file1,"%f\t",w[FIX(i,j,k)]);
        fprintf(file1,"%f\t",v[FIX(i,j-1,k)]);
        fprintf(file1,"%f\t",v[FIX(i,j,k)]);
        fprintf(file1,"\n");
        fprintf(file1,"%f\t",Dy*Dz*u[FIX(i-1,j,k)]);
        fprintf(file1,"%f\t",Dy*Dz*u[FIX(i,j,k)]);
        fprintf(file1,"%f\t",Dx*Dy*w[FIX(i,j,k-1)]);
        fprintf(file1,"%f\t",Dx*Dy*w[FIX(i,j,k)]);
        fprintf(file1,"%f\t",Dx*Dz*v[FIX(i,j-1,k)]);
        fprintf(file1,"%f\t",Dx*Dz*v[FIX(i,j,k)]);
        fprintf(file1,"\n");
        fprintf(file1,"\n");
    }

   }
  fclose(file1);
}

