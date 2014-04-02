#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "data_structure.h"
#include "solver_gs.h"
#include "boundary.h"
#include "utility.h"

FILE *file1;

/******************************************************************************
| Gauss-Seidel Solver
******************************************************************************/ 
void GS_P(PARA_DATA *para, REAL **var, int Type, REAL *x)
{
  REAL *as = var[AS], *aw = var[AW], *ae = var[AE], *an = var[AN];
  REAL *ap = var[AP], *af = var[AF], *ab = var[AB], *b = var[B];
  int imax = para->geom->imax, jmax= para->geom->jmax;
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);  
  int i, j, k, it=0;
  float tmp1,tmp2,residual=1;
  REAL *flagp = var[FLAGP],*flagu = var[FLAGU],*flagv = var[FLAGV],*flagw = var[FLAGW];


 
    while(residual>0.001 && it <5)

 // for(it=0; it<50; it++)
  {
         tmp1=0;
     tmp2=0.0000000001;
      it +=1;
     residual=0;



  for(i=1; i<=imax; i++)
      for(j=1; j<=jmax; j++)
        for(k=1; k<=kmax; k++)
    {
      if (flagp[FIX(i,j,k)]==1) continue;

          x[FIX(i,j,k)] = (  ae[FIX(i,j,k)]*x[FIX(i+1,j,k)] 
                          + aw[FIX(i,j,k)]*x[FIX(i-1,j,k)]
                          + an[FIX(i,j,k)]*x[FIX(i,j+1,k)]
                          + as[FIX(i,j,k)]*x[FIX(i,j-1,k)]
                          + af[FIX(i,j,k)]*x[FIX(i,j,k+1)]
                          + ab[FIX(i,j,k)]*x[FIX(i,j,k-1)]
                          + b[FIX(i,j,k)] ) / ap[FIX(i,j,k)];
    }

  for(i=imax; i>=1; i--)
      for(j=1; j<=jmax; j++)
        for(k=1; k<=kmax; k++)
    {
      if (flagp[FIX(i,j,k)]==1) continue;

          x[FIX(i,j,k)] = (  ae[FIX(i,j,k)]*x[FIX(i+1,j,k)] 
                          + aw[FIX(i,j,k)]*x[FIX(i-1,j,k)]
                          + an[FIX(i,j,k)]*x[FIX(i,j+1,k)]
                          + as[FIX(i,j,k)]*x[FIX(i,j-1,k)]
                          + af[FIX(i,j,k)]*x[FIX(i,j,k+1)]
                          + ab[FIX(i,j,k)]*x[FIX(i,j,k-1)]
                          + b[FIX(i,j,k)] ) / ap[FIX(i,j,k)];
    }

      for(j=1; j<=jmax; j++)
         for(i=1; i<=imax; i++)
            for(k=1; k<=kmax; k++)
    {
      if (flagp[FIX(i,j,k)]==1) continue;

          x[FIX(i,j,k)] = (  ae[FIX(i,j,k)]*x[FIX(i+1,j,k)] 
                          + aw[FIX(i,j,k)]*x[FIX(i-1,j,k)]
                          + an[FIX(i,j,k)]*x[FIX(i,j+1,k)]
                          + as[FIX(i,j,k)]*x[FIX(i,j-1,k)]
                          + af[FIX(i,j,k)]*x[FIX(i,j,k+1)]
                          + ab[FIX(i,j,k)]*x[FIX(i,j,k-1)]
                          + b[FIX(i,j,k)] ) / ap[FIX(i,j,k)];
    }
      for(j=jmax; j>=1; j--)
         for(i=1; i<=imax; i++)
            for(k=1; k<=kmax; k++)
    {
      if (flagp[FIX(i,j,k)]==1) continue;

          x[FIX(i,j,k)] = (  ae[FIX(i,j,k)]*x[FIX(i+1,j,k)] 
                          + aw[FIX(i,j,k)]*x[FIX(i-1,j,k)]
                          + an[FIX(i,j,k)]*x[FIX(i,j+1,k)]
                          + as[FIX(i,j,k)]*x[FIX(i,j-1,k)]
                          + af[FIX(i,j,k)]*x[FIX(i,j,k+1)]
                          + ab[FIX(i,j,k)]*x[FIX(i,j,k-1)]
                          + b[FIX(i,j,k)] ) / ap[FIX(i,j,k)];
    }
 

  //for(k=1; k<=kmax; k++)
  //   for(i=1; i<=imax; i++)
  //     for(j=1; j<=jmax; j++)
  //  {
  //    if (flagp[FIX(i,j,k)]==1) continue;

  //        x[FIX(i,j,k)] = (  ae[FIX(i,j,k)]*x[FIX(i+1,j,k)] 
  //                        + aw[FIX(i,j,k)]*x[FIX(i-1,j,k)]
  //                        + an[FIX(i,j,k)]*x[FIX(i,j+1,k)]
  //                        + as[FIX(i,j,k)]*x[FIX(i,j-1,k)]
  //                        + af[FIX(i,j,k)]*x[FIX(i,j,k+1)]
  //                        + ab[FIX(i,j,k)]*x[FIX(i,j,k-1)]
  //                        + b[FIX(i,j,k)] ) / ap[FIX(i,j,k)];
  //  }
  // for(k=kmax; k>=1; k--)
  //   for(i=1; i<=imax; i++)
  //     for(j=1; j<=jmax; j++)
  //  {
  //    if (flagp[FIX(i,j,k)]==1) continue;

  //        x[FIX(i,j,k)] = (  ae[FIX(i,j,k)]*x[FIX(i+1,j,k)] 
  //                        + aw[FIX(i,j,k)]*x[FIX(i-1,j,k)]
  //                        + an[FIX(i,j,k)]*x[FIX(i,j+1,k)]
  //                        + as[FIX(i,j,k)]*x[FIX(i,j-1,k)]
  //                        + af[FIX(i,j,k)]*x[FIX(i,j,k+1)]
  //                        + ab[FIX(i,j,k)]*x[FIX(i,j,k-1)]
  //                        + b[FIX(i,j,k)] ) / ap[FIX(i,j,k)];
  //  }
 
  //for(i=1; i<=imax; i++)
  //    for(j=jmax; j>=1; j--)
  //      for(k=kmax; k>=1; k--)
  //  {
  //    if (flagp[FIX(i,j,k)]==1) continue;

  //        x[FIX(i,j,k)] = (  ae[FIX(i,j,k)]*x[FIX(i+1,j,k)] 
  //                        + aw[FIX(i,j,k)]*x[FIX(i-1,j,k)]
  //                        + an[FIX(i,j,k)]*x[FIX(i,j+1,k)]
  //                        + as[FIX(i,j,k)]*x[FIX(i,j-1,k)]
  //                        + af[FIX(i,j,k)]*x[FIX(i,j,k+1)]
  //                        + ab[FIX(i,j,k)]*x[FIX(i,j,k-1)]
  //                        + b[FIX(i,j,k)] ) / ap[FIX(i,j,k)];
  //  }

  //for(i=imax; i>=1; i--)
  //    for(j=jmax; j>=1; j--)
  //      for(k=kmax; k>=1; k--)
  //  {
  //    if (flagp[FIX(i,j,k)]==1) continue;

  //        x[FIX(i,j,k)] = (  ae[FIX(i,j,k)]*x[FIX(i+1,j,k)] 
  //                        + aw[FIX(i,j,k)]*x[FIX(i-1,j,k)]
  //                        + an[FIX(i,j,k)]*x[FIX(i,j+1,k)]
  //                        + as[FIX(i,j,k)]*x[FIX(i,j-1,k)]
  //                        + af[FIX(i,j,k)]*x[FIX(i,j,k+1)]
  //                        + ab[FIX(i,j,k)]*x[FIX(i,j,k-1)]
  //                        + b[FIX(i,j,k)] ) / ap[FIX(i,j,k)];
  //  }

  //    for(j=1; j<=jmax; j++)
  //       for(i=imax; i>=1; i--)
  //         for(k=kmax; k>=1; k--)
  //  {
  //    if (flagp[FIX(i,j,k)]==1) continue;

  //        x[FIX(i,j,k)] = (  ae[FIX(i,j,k)]*x[FIX(i+1,j,k)] 
  //                        + aw[FIX(i,j,k)]*x[FIX(i-1,j,k)]
  //                        + an[FIX(i,j,k)]*x[FIX(i,j+1,k)]
  //                        + as[FIX(i,j,k)]*x[FIX(i,j-1,k)]
  //                        + af[FIX(i,j,k)]*x[FIX(i,j,k+1)]
  //                        + ab[FIX(i,j,k)]*x[FIX(i,j,k-1)]
  //                        + b[FIX(i,j,k)] ) / ap[FIX(i,j,k)];
  //  }
  //    for(j=jmax; j>=1; j--)
  //       for(i=imax; i>=1; i--)
  //         for(k=kmax; k>=1; k--)
  //  {
  //    if (flagp[FIX(i,j,k)]==1) continue;

  //        x[FIX(i,j,k)] = (  ae[FIX(i,j,k)]*x[FIX(i+1,j,k)] 
  //                        + aw[FIX(i,j,k)]*x[FIX(i-1,j,k)]
  //                        + an[FIX(i,j,k)]*x[FIX(i,j+1,k)]
  //                        + as[FIX(i,j,k)]*x[FIX(i,j-1,k)]
  //                        + af[FIX(i,j,k)]*x[FIX(i,j,k+1)]
  //                        + ab[FIX(i,j,k)]*x[FIX(i,j,k-1)]
  //                        + b[FIX(i,j,k)] ) / ap[FIX(i,j,k)];
  //  }
 

  //for(k=1; k<=kmax; k++)
  //       for(i=imax; i>=1; i--)
  //         for(j=jmax; j>=1; j--)
  //  {
  //    if (flagp[FIX(i,j,k)]==1) continue;

  //        x[FIX(i,j,k)] = (  ae[FIX(i,j,k)]*x[FIX(i+1,j,k)] 
  //                        + aw[FIX(i,j,k)]*x[FIX(i-1,j,k)]
  //                        + an[FIX(i,j,k)]*x[FIX(i,j+1,k)]
  //                        + as[FIX(i,j,k)]*x[FIX(i,j-1,k)]
  //                        + af[FIX(i,j,k)]*x[FIX(i,j,k+1)]
  //                        + ab[FIX(i,j,k)]*x[FIX(i,j,k-1)]
  //                        + b[FIX(i,j,k)] ) / ap[FIX(i,j,k)];
  //  }
  // for(k=kmax; k>=1; k--)
  //       for(i=imax; i>=1; i--)
  //         for(j=jmax; j>=1; j--)
  //  {
  //    if (flagp[FIX(i,j,k)]==1) continue;

  //        x[FIX(i,j,k)] = (  ae[FIX(i,j,k)]*x[FIX(i+1,j,k)] 
  //                        + aw[FIX(i,j,k)]*x[FIX(i-1,j,k)]
  //                        + an[FIX(i,j,k)]*x[FIX(i,j+1,k)]
  //                        + as[FIX(i,j,k)]*x[FIX(i,j-1,k)]
  //                        + af[FIX(i,j,k)]*x[FIX(i,j,k+1)]
  //                        + ab[FIX(i,j,k)]*x[FIX(i,j,k-1)]
  //                        + b[FIX(i,j,k)] ) / ap[FIX(i,j,k)];
  //  }

   FOR_EACH_CELL
          if (flagp[FIX(i,j,k)]==1) continue;
    tmp1 += (float) fabs(ap[FIX(i,j,k)]*x[FIX(i,j,k)] 
          - ae[FIX(i,j,k)]*x[FIX(i+1,j,k)] - aw[FIX(i,j,k)]*x[FIX(i-1,j,k)]
          - an[FIX(i,j,k)]*x[FIX(i,j+1,k)] - as[FIX(i,j,k)]*x[FIX(i,j-1,k)]
          - af[FIX(i,j,k)]*x[FIX(i,j,k+1)] - ab[FIX(i,j,k)]*x[FIX(i,j,k-1)]
          - b[FIX(i,j,k)]);
    tmp2 += (float) fabs(ap[FIX(i,j,k)]*x[FIX(i,j,k)]);


  //residual= residual> (tmp1/tmp2)?residual:(tmp1/tmp2);

   END_FOR
     residual  = (tmp1 /tmp2);
   
  

  //printf("it=%d\n",it);

  //fprintf (file1,"%d\t%.3g\n",it,residual);

      para->prob->resp=residual;

 // printf ("residual=%.3g\n",residual);

      para->prob->iterp +=1;      

   
}

 //  fclose(file1);

} // End of Gauss-Seidel( )

void Gauss_Seidel(PARA_DATA *para, REAL **var, REAL *flag, REAL *x)
{
  REAL *as = var[AS], *aw = var[AW], *ae = var[AE], *an = var[AN];
  REAL *ap = var[AP], *af = var[AF], *ab = var[AB], *b = var[B];
  int imax = para->geom->imax, jmax= para->geom->jmax;
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);  
  int i, j, k, it=0;
  float tmp1,tmp2,residual=1;

 while(residual>0.001 && it<5) //residual>0.001 &&
  // for(it=0; it<5; it++)
  {

         tmp1=0;
     tmp2=0.0000000001;
     it +=1;
     residual=0;

  for(i=1; i<=imax; i++)
      for(j=1; j<=jmax; j++)
        for(k=1; k<=kmax; k++)
    {
      if (flag[FIX(i,j,k)]>=0) continue;

          x[FIX(i,j,k)] = (  ae[FIX(i,j,k)]*x[FIX(i+1,j,k)] 
                          + aw[FIX(i,j,k)]*x[FIX(i-1,j,k)]
                          + an[FIX(i,j,k)]*x[FIX(i,j+1,k)]
                          + as[FIX(i,j,k)]*x[FIX(i,j-1,k)]
                          + af[FIX(i,j,k)]*x[FIX(i,j,k+1)]
                          + ab[FIX(i,j,k)]*x[FIX(i,j,k-1)]
                          + b[FIX(i,j,k)] ) / ap[FIX(i,j,k)];
    }

  for(i=imax; i>=1; i--)
      for(j=1; j<=jmax; j++)
        for(k=1; k<=kmax; k++)
    {
      if (flag[FIX(i,j,k)]>=0) continue;

          x[FIX(i,j,k)] = (  ae[FIX(i,j,k)]*x[FIX(i+1,j,k)] 
                          + aw[FIX(i,j,k)]*x[FIX(i-1,j,k)]
                          + an[FIX(i,j,k)]*x[FIX(i,j+1,k)]
                          + as[FIX(i,j,k)]*x[FIX(i,j-1,k)]
                          + af[FIX(i,j,k)]*x[FIX(i,j,k+1)]
                          + ab[FIX(i,j,k)]*x[FIX(i,j,k-1)]
                          + b[FIX(i,j,k)] ) / ap[FIX(i,j,k)];
    }
      for(j=1; j<=jmax; j++)
         for(i=1; i<=imax; i++)
            for(k=1; k<=kmax; k++)
    {
      if (flag[FIX(i,j,k)]>=0) continue;

          x[FIX(i,j,k)] = (  ae[FIX(i,j,k)]*x[FIX(i+1,j,k)] 
                          + aw[FIX(i,j,k)]*x[FIX(i-1,j,k)]
                          + an[FIX(i,j,k)]*x[FIX(i,j+1,k)]
                          + as[FIX(i,j,k)]*x[FIX(i,j-1,k)]
                          + af[FIX(i,j,k)]*x[FIX(i,j,k+1)]
                          + ab[FIX(i,j,k)]*x[FIX(i,j,k-1)]
                          + b[FIX(i,j,k)] ) / ap[FIX(i,j,k)];
    }
      for(j=jmax; j>=1; j--)
         for(i=1; i<=imax; i++)
            for(k=1; k<=kmax; k++)
    {
      if (flag[FIX(i,j,k)]>=0) continue;

          x[FIX(i,j,k)] = (  ae[FIX(i,j,k)]*x[FIX(i+1,j,k)] 
                          + aw[FIX(i,j,k)]*x[FIX(i-1,j,k)]
                          + an[FIX(i,j,k)]*x[FIX(i,j+1,k)]
                          + as[FIX(i,j,k)]*x[FIX(i,j-1,k)]
                          + af[FIX(i,j,k)]*x[FIX(i,j,k+1)]
                          + ab[FIX(i,j,k)]*x[FIX(i,j,k-1)]
                          + b[FIX(i,j,k)] ) / ap[FIX(i,j,k)];
    }
 

  //for(k=1; k<=kmax; k++)
  //   for(i=1; i<=imax; i++)
  //     for(j=1; j<=jmax; j++)
  //  {
  //    if (flag[FIX(i,j,k)]>=0) continue;

  //        x[FIX(i,j,k)] = (  ae[FIX(i,j,k)]*x[FIX(i+1,j,k)] 
  //                        + aw[FIX(i,j,k)]*x[FIX(i-1,j,k)]
  //                        + an[FIX(i,j,k)]*x[FIX(i,j+1,k)]
  //                        + as[FIX(i,j,k)]*x[FIX(i,j-1,k)]
  //                        + af[FIX(i,j,k)]*x[FIX(i,j,k+1)]
  //                        + ab[FIX(i,j,k)]*x[FIX(i,j,k-1)]
  //                        + b[FIX(i,j,k)] ) / ap[FIX(i,j,k)];
  //  }
  // for(k=kmax; k>=1; k--)
  //   for(i=1; i<=imax; i++)
  //     for(j=1; j<=jmax; j++)
  //  {
  //    if (flag[FIX(i,j,k)]>=0) continue;

  //        x[FIX(i,j,k)] = (  ae[FIX(i,j,k)]*x[FIX(i+1,j,k)] 
  //                        + aw[FIX(i,j,k)]*x[FIX(i-1,j,k)]
  //                        + an[FIX(i,j,k)]*x[FIX(i,j+1,k)]
  //                        + as[FIX(i,j,k)]*x[FIX(i,j-1,k)]
  //                        + af[FIX(i,j,k)]*x[FIX(i,j,k+1)]
  //                        + ab[FIX(i,j,k)]*x[FIX(i,j,k-1)]
  //                        + b[FIX(i,j,k)] ) / ap[FIX(i,j,k)];
  //  }
  //for(i=1; i<=imax; i++)
  //    for(j=jmax; j>=1; j--)
  //      for(k=kmax; k>=1; k--)
  //  {
  //    if (flag[FIX(i,j,k)]>=0) continue;

  //        x[FIX(i,j,k)] = (  ae[FIX(i,j,k)]*x[FIX(i+1,j,k)] 
  //                        + aw[FIX(i,j,k)]*x[FIX(i-1,j,k)]
  //                        + an[FIX(i,j,k)]*x[FIX(i,j+1,k)]
  //                        + as[FIX(i,j,k)]*x[FIX(i,j-1,k)]
  //                        + af[FIX(i,j,k)]*x[FIX(i,j,k+1)]
  //                        + ab[FIX(i,j,k)]*x[FIX(i,j,k-1)]
  //                        + b[FIX(i,j,k)] ) / ap[FIX(i,j,k)];
  //  }

  //for(i=imax; i>=1; i--)
  //    for(j=jmax; j>=1; j--)
  //      for(k=kmax; k>=1; k--)
  //  {
  //    if (flag[FIX(i,j,k)]>=0) continue;

  //        x[FIX(i,j,k)] = (  ae[FIX(i,j,k)]*x[FIX(i+1,j,k)] 
  //                        + aw[FIX(i,j,k)]*x[FIX(i-1,j,k)]
  //                        + an[FIX(i,j,k)]*x[FIX(i,j+1,k)]
  //                        + as[FIX(i,j,k)]*x[FIX(i,j-1,k)]
  //                        + af[FIX(i,j,k)]*x[FIX(i,j,k+1)]
  //                        + ab[FIX(i,j,k)]*x[FIX(i,j,k-1)]
  //                        + b[FIX(i,j,k)] ) / ap[FIX(i,j,k)];
  //  }

  //    for(j=1; j<=jmax; j++)
  //       for(i=imax; i>=1; i--)
  //         for(k=kmax; k>=1; k--)
  //  {
  //    if (flag[FIX(i,j,k)]>=0) continue;

  //        x[FIX(i,j,k)] = (  ae[FIX(i,j,k)]*x[FIX(i+1,j,k)] 
  //                        + aw[FIX(i,j,k)]*x[FIX(i-1,j,k)]
  //                        + an[FIX(i,j,k)]*x[FIX(i,j+1,k)]
  //                        + as[FIX(i,j,k)]*x[FIX(i,j-1,k)]
  //                        + af[FIX(i,j,k)]*x[FIX(i,j,k+1)]
  //                        + ab[FIX(i,j,k)]*x[FIX(i,j,k-1)]
  //                        + b[FIX(i,j,k)] ) / ap[FIX(i,j,k)];
  //  }
  //    for(j=jmax; j>=1; j--)
  //       for(i=imax; i>=1; i--)
  //         for(k=kmax; k>=1; k--)
  //  {
  //    if (flag[FIX(i,j,k)]>=0) continue;

  //        x[FIX(i,j,k)] = (  ae[FIX(i,j,k)]*x[FIX(i+1,j,k)] 
  //                        + aw[FIX(i,j,k)]*x[FIX(i-1,j,k)]
  //                        + an[FIX(i,j,k)]*x[FIX(i,j+1,k)]
  //                        + as[FIX(i,j,k)]*x[FIX(i,j-1,k)]
  //                        + af[FIX(i,j,k)]*x[FIX(i,j,k+1)]
  //                        + ab[FIX(i,j,k)]*x[FIX(i,j,k-1)]
  //                        + b[FIX(i,j,k)] ) / ap[FIX(i,j,k)];
  //  }
 

  //for(k=1; k<=kmax; k++)
  //       for(i=imax; i>=1; i--)
  //         for(j=jmax; j>=1; j--)
  //  {
  //    if (flag[FIX(i,j,k)]>=0) continue;

  //        x[FIX(i,j,k)] = (  ae[FIX(i,j,k)]*x[FIX(i+1,j,k)] 
  //                        + aw[FIX(i,j,k)]*x[FIX(i-1,j,k)]
  //                        + an[FIX(i,j,k)]*x[FIX(i,j+1,k)]
  //                        + as[FIX(i,j,k)]*x[FIX(i,j-1,k)]
  //                        + af[FIX(i,j,k)]*x[FIX(i,j,k+1)]
  //                        + ab[FIX(i,j,k)]*x[FIX(i,j,k-1)]
  //                        + b[FIX(i,j,k)] ) / ap[FIX(i,j,k)];
  //  }
  // for(k=kmax; k>=1; k--)
  //       for(i=imax; i>=1; i--)
  //         for(j=jmax; j>=1; j--)
  //  {
  //    if (flag[FIX(i,j,k)]>=0) continue;

  //        x[FIX(i,j,k)] = (  ae[FIX(i,j,k)]*x[FIX(i+1,j,k)] 
  //                        + aw[FIX(i,j,k)]*x[FIX(i-1,j,k)]
  //                        + an[FIX(i,j,k)]*x[FIX(i,j+1,k)]
  //                        + as[FIX(i,j,k)]*x[FIX(i,j-1,k)]
  //                        + af[FIX(i,j,k)]*x[FIX(i,j,k+1)]
  //                        + ab[FIX(i,j,k)]*x[FIX(i,j,k-1)]
  //                        + b[FIX(i,j,k)] ) / ap[FIX(i,j,k)];
  //  } 
 
     
   FOR_EACH_CELL
                  if (flag[FIX(i,j,k)]>=0) continue;
    tmp1 += (float) fabs(ap[FIX(i,j,k)]*x[FIX(i,j,k)] 
          - ae[FIX(i,j,k)]*x[FIX(i+1,j,k)] - aw[FIX(i,j,k)]*x[FIX(i-1,j,k)]
          - an[FIX(i,j,k)]*x[FIX(i,j+1,k)] - as[FIX(i,j,k)]*x[FIX(i,j-1,k)]
          - af[FIX(i,j,k)]*x[FIX(i,j,k+1)] - ab[FIX(i,j,k)]*x[FIX(i,j,k-1)]
          - b[FIX(i,j,k)]);
    tmp2 += (float) fabs(ap[FIX(i,j,k)]*x[FIX(i,j,k)]);


  //residual= residual> (tmp1/tmp2)?residual:(tmp1/tmp2);

   END_FOR
     residual  = (tmp1 /tmp2);
   

   para->prob->iteru +=1;



}


}// End of Gauss-Seidel( )



void Gauss_Seidel_simple(PARA_DATA *para, REAL **var, REAL *flag, REAL *x)
{
  REAL *as = var[AS], *aw = var[AW], *ae = var[AE], *an = var[AN];
  REAL *ap = var[AP], *af = var[AF], *ab = var[AB], *b = var[B];
  int imax = para->geom->imax, jmax= para->geom->jmax;
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);  
  int i, j, k, it=0;

  for(it=0; it<5; it++)
  {

    for(i=1; i<=imax; i++)
      for(j=1; j<=jmax; j++)
        for(k=1; k<=kmax; k++)
    {
                  if (flag[FIX(i,j,k)]==1) continue;
          x[FIX(i,j,k)] = (  ae[FIX(i,j,k)]*x[FIX(i+1,j,k)] 
                          + aw[FIX(i,j,k)]*x[FIX(i-1,j,k)]
                          + an[FIX(i,j,k)]*x[FIX(i,j+1,k)]
                          + as[FIX(i,j,k)]*x[FIX(i,j-1,k)]
                          + af[FIX(i,j,k)]*x[FIX(i,j,k+1)]
                          + ab[FIX(i,j,k)]*x[FIX(i,j,k-1)]
                          + b[FIX(i,j,k)] ) / ap[FIX(i,j,k)];
    }
    
    for(i=imax;i>=1;i--)
      for(j=jmax;j>=1;j--)
        for(k=kmax;k>=1;k--)
    {
                 if (flag[FIX(i,j,k)]==1) continue;
          x[FIX(i,j,k)] = (  ae[FIX(i,j,k)]*x[FIX(i+1,j,k)] 
                          + aw[FIX(i,j,k)]*x[FIX(i-1,j,k)]
                          + an[FIX(i,j,k)]*x[FIX(i,j+1,k)]
                          + as[FIX(i,j,k)]*x[FIX(i,j-1,k)]
                          + af[FIX(i,j,k)]*x[FIX(i,j,k+1)]
                          + ab[FIX(i,j,k)]*x[FIX(i,j,k-1)]
                          + b[FIX(i,j,k)] ) / ap[FIX(i,j,k)]; 
    }

  }
   


}// End of Gauss-Seidel( )

