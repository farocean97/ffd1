#define _CRT_SECURE_NO_DEPRECATE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "data_structure.h"
#include "read_data.h"


FILE *file_params;

/******************************************************************************
| Read the data from input_displacement_li_coarse .cfd
******************************************************************************/

int read_max(PARA_DATA *para, REAL **var) {  
  int imax,jmax,kmax;
  float Lx,Ly,Lz;
  char  string[400];

  if( (file_params=fopen("input_displacement_li_coarse.cfd","r")) == NULL ) {
    fprintf(stderr,"Error:can not open error file!\n");
    return 0;
  }
  fgets(string, 400, file_params);
  sscanf(string,"%f%f%f",&Lx,&Ly,&Lz);

  fgets(string, 400, file_params);
  sscanf(string,"%d%d%d", &imax,&jmax,&kmax);
  fclose(file_params);

  para->geom->imax=imax;
  para->geom->jmax=jmax;
  para->geom->kmax=kmax;
  printf("%d\t%d\t%d\n",imax,jmax,kmax);

  para->geom->Lx=Lx;
  para->geom->Ly=Ly;
  para->geom->Lz=Lz;

  return 1;

}

int read_input(PARA_DATA *para, REAL **var, int **BINDEX) {
  int i,j, k;
  REAL tempx=0.0f, tempy=0.0f, tempz=0.0f;
  REAL Lx = para->geom->Lx;
  REAL Ly = para->geom->Ly;
  REAL Lz = para->geom->Lz;
  REAL *gx = var[GX], *gy = var[GY], *gz = var[GZ];
  REAL *x = var[X], *y = var[Y], *z = var[Z];
  int IWWALL,IEWALL,ISWALL,INWALL,IBWALL,ITWALL;
  int NBIN, NBOUT, NBL, NW;
  int NBUS,NBS;
  int SI,SJ,SK,EI,EJ,EK,NI,NJ,NK,FLTMP;
  float TMP,MASS,U,V,W;
  int restart;
  float density,nu,cp,gravx,gravy,gravz,beta,trefmax,spec;
  float t_start,t_delta,t_total;
  int imax = para->geom->imax;
  int jmax = para->geom->jmax;
  int kmax = para->geom->kmax;
  int index=1;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2); 
  int zone_num=0;
  char *temp, string[400];
  float *delx,*dely,*delz;
  REAL *flagp = var[FLAGP],*flagu = var[FLAGU];
  REAL *flagv = var[FLAGV],*flagw = var[FLAGW];

  if( (file_params=fopen("input_displacement_li_coarse.cfd","r")) == NULL) {
    fprintf(stderr,"Error:can not open error file!\n");
    return 0;
  }
    
  delx = (float *) malloc ((imax+2)*sizeof(float));
  dely = (float *) malloc ((jmax+2)*sizeof(float));
  delz = (float *) malloc ((kmax+2)*sizeof(float));

  if ( !delx || !dely ||!delz ) {
    fprintf ( stderr, "cannot allocate data\n" );
    return ( 0 );
  }

  delx[0]=0;
  dely[0]=0;
  delz[0]=0;

  temp = fgets(string, 400, file_params);
  temp = fgets(string, 400, file_params);

  for(i=1;i<=imax;i++)
  fscanf(file_params,"%f" ,&delx[i]); 
  fscanf(file_params,"\n");

  for(j=1;j<=jmax;j++)
  fscanf(file_params,"%f" ,&dely[j]); 
  fscanf(file_params,"\n");

  for(k=1;k<=kmax;k++)
  fscanf(file_params,"%f" ,&delz[k]); 
  fscanf(file_params,"\n");

  for(i=0;i<=imax+1;i++) {
    tempx += delx[i];
    if(i==imax) tempx=Lx;
    if(i>imax)  tempx=Lx+0.000001f;
    for(j=0;j<=jmax+1;j++)
      for(k=0;k<=kmax+1;k++) {
        var[GX][FIX(i,j,k)]=tempx;
      }
  }

  for(j=0;j<=jmax+1;j++)  {
    tempy += dely[j];
    if(j==jmax) tempy=Ly;
    if(j>jmax)  tempy=Ly+0.000001f;
      for(i=0;i<=imax+1;i++)
        for(k=0;k<=kmax+1;k++) {
          var[GY][FIX(i,j,k)]=tempy;
        }
  }

  for(k=0;k<=kmax+1;k++) {
    tempz += delz[k];
    if(k==kmax) tempz=Lz;
    if(k>kmax)  tempz=Lz+0.000001f;
    for(i=0;i<=imax+1;i++)
      for(j=0;j<=jmax+1;j++) {
        var[GZ][FIX(i,j,k)]=tempz;
      }
  }

  FOR_ALL_CELL

    if(i<1)  x[FIX(i,j,k)]= 0;
    else if(i>imax) x[FIX(i,j,k)]= Lx;
    else  x[FIX(i,j,k)]= 0.5f* (gx[FIX(i,j,k)]+gx[FIX(i-1,j,k)]);
   
    if(j<1)  y[FIX(i,j,k)]= 0;
    else if(j>jmax) y[FIX(i,j,k)]= Ly;
    else y[FIX(i,j,k)]= 0.5f* (gy[FIX(i,j,k)]+gy[FIX(i,j-1,k)]);

    if(k<1)  z[FIX(i,j,k)]= 0;
    else if(k>kmax) z[FIX(i,j,k)]= Lz;
    else z[FIX(i,j,k)]= 0.5f* (gz[FIX(i,j,k)]+gz[FIX(i,j,k-1)]);

  END_FOR

/**********************************      INLET      ************************************************/
  fgets(string, 400, file_params);
  sscanf(string,"%d%d%d%d%d%d",&IWWALL,&IEWALL,&ISWALL,&INWALL,&IBWALL,&ITWALL); 
  fgets(string, 400, file_params);
  sscanf(string,"%d",&NBIN); 

  if(NBIN != 0) {
    for(i=1;i<=NBIN;i++) {
      fgets(string, 400, file_params);
      sscanf(string,"%d%d%d%d%d%d%f%f%f%f%f",&SI,&SJ,&SK ,&NI,&NJ,
                                              &NK ,&TMP ,&MASS ,&U ,&V ,&W );
      if(SI !=0) SI -= 1;
      if(SJ !=0) SJ -= 1;
      if(SK !=0) SK -= 1;
      EI=SI+NI;
      EJ=SJ+NJ;
      EK=SK+NK;
  
      BINDEX[7][zone_num+i]= SI;
      BINDEX[8][zone_num+i]= EI;
      BINDEX[9][zone_num+i]= SJ;
      BINDEX[10][zone_num+i]= EJ;
      BINDEX[11][zone_num+i]= SK;
      BINDEX[12][zone_num+i]= EK;
      para->bc->u_bc[zone_num+i]=U;
      para->bc->v_bc[zone_num+i]=V;
      para->bc->w_bc[zone_num+i]=W;
      para->bc->t_bc[zone_num+i]=TMP;
      para->bc->um_bc[zone_num+i]=MASS;
    }
  }//end NBIN
      
  zone_num += NBIN;
  para->geom->zone_inlet=zone_num;

   
/**********************************      OUTLET      ************************************************/
   
  fgets(string, 400, file_params);
  sscanf(string,"%d",&NBOUT); 
  para->bc->NBOUT=NBOUT;

  if(NBOUT !=0) {
    for(i=1;i<=NBOUT;i++)  {
      fgets(string, 400, file_params);
      sscanf(string,"%d%d%d%d%d%d%f%f%f%f%f",&SI,&SJ,&SK ,&NI,&NJ,
                                             &NK ,&TMP ,&MASS ,&U ,&V ,&W );
      if(SI !=0) SI -= 1;
      if(SJ !=0) SJ -= 1;
      if(SK !=0) SK -= 1;
      EI=SI+NI;
      EJ=SJ+NJ;
      EK=SK+NK;

      BINDEX[7][zone_num+i]= SI;
      BINDEX[8][zone_num+i]= EI;
      BINDEX[9][zone_num+i]= SJ;
      BINDEX[10][zone_num+i]= EJ;
      BINDEX[11][zone_num+i]= SK;
      BINDEX[12][zone_num+i]= EK;
      para->bc->u_bc[zone_num+i]=U;
      para->bc->v_bc[zone_num+i]=V;
      para->bc->w_bc[zone_num+i]=W;
      para->bc->t_bc[zone_num+i]=TMP;
    }
  }

  zone_num +=NBOUT;
  para->geom->zone_outlet=zone_num;


/**********************************      BLOCKAGE     ************************************************/

  fgets(string, 400, file_params);
  sscanf(string,"%d",&NBL); 

  if(NBL !=0) {
    for(i=1;i<=NBL;i++)  {
      fgets(string, 400, file_params);
      sscanf(string,"%d%d%d%d%d%d%d%f",&SI,&SJ,&SK ,&NI,
                                       &NJ ,&NK ,&FLTMP, &TMP);
      if(SI !=0) SI -= 1;
      if(SJ !=0) SJ -= 1;
      if(SK !=0) SK -= 1;
      EI=SI+NI;
      EJ=SJ+NJ;
      EK=SK+NK;

      BINDEX[7][zone_num+i]= SI;
      BINDEX[8][zone_num+i]= EI;
      BINDEX[9][zone_num+i]= SJ;
      BINDEX[10][zone_num+i]= EJ;
      BINDEX[11][zone_num+i]= SK;
      BINDEX[12][zone_num+i]= EK;
      para->bc->t_bc[zone_num+i]=TMP;
      para->bc->fltmp[zone_num+i]=FLTMP; //FLTMP=1 Const temp, FLTMP=0, Const heatflux
    }
  }
  zone_num += NBL;
  para->geom->zone_bl=zone_num;

/**********************************      Wall  ************************************************/

  fgets(string, 400, file_params);
  sscanf(string,"%d",&NW); 

  if(NW !=0) {
    for(i=1;i<=NW;i++)  {
      fgets(string, 400, file_params);
      sscanf(string,"%d%d%d%d%d%d%d%f",&SI,&SJ,&SK ,&NI,&NJ,
                                       &NK ,&FLTMP, &TMP);
      if(SI !=0) SI -= 1;
      if(SJ !=0) SJ -= 1;
      if(SK !=0) SK -= 1;
      EI=SI+NI;
      EJ=SJ+NJ;
      EK=SK+NK;

      BINDEX[7][zone_num+i]= SI;
      BINDEX[8][zone_num+i]= EI;
      BINDEX[9][zone_num+i]= SJ;
      BINDEX[10][zone_num+i]= EJ;
      BINDEX[11][zone_num+i]= SK;
      BINDEX[12][zone_num+i]= EK;
      para->bc->t_bc[zone_num+i]=TMP;
      para->bc->fltmp[zone_num+i]=FLTMP;
    }
  }
  zone_num += NW;
  para->geom->zone_wall=zone_num;

   /**********************************      Source     ************************************************/

  fgets(string, 400, file_params);
  sscanf(string,"%d",&NBUS); 

  if(NBUS !=0) {
    for(i=1;i<=NBUS;i++) {
      fgets(string, 400, file_params);
      sscanf(string,"%d%d%d%d%d%d%f",&SI,&SJ,&SK ,&NI,&NJ ,&NK , &TMP);

      if(SI !=0) SI -= 1;
      if(SJ !=0) SJ -= 1;
      if(SK !=0) SK -= 1;
      EI=SI+NI;
      EJ=SJ+NJ;
      EK=SK+NK;

      BINDEX[7][zone_num+i]= SI;
      BINDEX[8][zone_num+i]= EI;
      BINDEX[9][zone_num+i]= SJ;
      BINDEX[10][zone_num+i]= EJ;
      BINDEX[11][zone_num+i]= SK;
      BINDEX[12][zone_num+i]= EK;
      para->bc->t_bc[zone_num+i]=TMP/(REAL)(NI*NJ*NK);
      para->bc->qs +=TMP;
    }
  }
  zone_num += NBUS;
  para->geom->zone_us=zone_num;
  //para->geom->zone_num = zone_num;

  /**********************************      Openings in the walls     ************************************************/

  fgets(string, 400, file_params);
  sscanf(string,"%d",&NBS); 

  if(NBS !=0) {
    for(i=1;i<=NBS;i++) {
      fgets(string, 400, file_params);
      sscanf(string,"%d%d%d%d%d%d",&SI,&SJ,&SK ,&NI,&NJ ,&NK);

      if(SI !=0) SI -= 1;
      if(SJ !=0) SJ -= 1;
      if(SK !=0) SK -= 1;
      EI=SI+NI;
      EJ=SJ+NJ;
      EK=SK+NK;

      BINDEX[7][zone_num+i]= SI;
      BINDEX[8][zone_num+i]= EI;
      BINDEX[9][zone_num+i]= SJ;
      BINDEX[10][zone_num+i]= EJ;
      BINDEX[11][zone_num+i]= SK;
      BINDEX[12][zone_num+i]= EK;
    }
  }
  zone_num += NBS;
  para->geom->zone_s=zone_num;
  para->geom->zone_num = zone_num;

 
/**********************************     Other parameter  ************************************************/

  temp = fgets(string, 400, file_params); //maximum iteration
  temp = fgets(string, 400, file_params); //convergence rate
  temp = fgets(string, 400, file_params); //Turbulence model
  temp = fgets(string, 400, file_params); //initial value
  temp = fgets(string, 400, file_params); //minimum value
  temp = fgets(string, 400, file_params); //maximum value
  temp = fgets(string, 400, file_params); //fts value
  temp = fgets(string, 400, file_params); //under relaxation
  temp = fgets(string, 400, file_params); //reference point
  temp = fgets(string, 400, file_params); //monitering point

  fgets(string, 400, file_params);
  sscanf(string,"%d",&restart);
  para->solv->read_file=restart;

  temp = fgets(string, 400, file_params); //print frequency
  temp = fgets(string, 400, file_params); //Pressure variable Y/N
  temp = fgets(string, 400, file_params); //Steady state, hall_partition_geom.

  fgets(string, 400, file_params);
  sscanf(string,"%f %f %f %f %f %f %f %f %f",&density,&nu,&cp,&gravx,
                                            &gravy,&gravz,&beta,&trefmax,&spec);
  para->prob->rho=density;
  para->prob->nu=nu;
  para->prob->cond=cp;
  para->prob->gravx=gravx;
  para->prob->gravy=gravy;
  para->prob->gravz=gravz;
  para->prob->beta=beta;
  para->prob->trefmax=trefmax;
  para->prob->spec=spec;

  fgets(string, 400, file_params);
  sscanf(string,"%f %f %f",&t_start,&t_delta,&t_total);
  para->mytime->t_start=t_start;
  para->mytime->dt=t_delta;
  para->mytime->t_output=t_total;

  temp = fgets(string, 400, file_params); //prandtl

  fclose(file_params);

  free(delx);
  free(dely);
  free(delz);
 
  return 1;

} // End of read_dara()


int read_zeroone(PARA_DATA *para, REAL **var, int **BINDEX) {
  int i,j, k;
  int delcount=0;
  int mark;
  int imax = para->geom->imax;
  int jmax = para->geom->jmax;
  int kmax = para->geom->kmax;
  int index = para->geom->index[1]+1;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  REAL *flagp = var[FLAGP],*flagu = var[FLAGU];
  REAL *flagv = var[FLAGV],*flagw = var[FLAGW];

  if( (file_params=fopen("zeroone_hall_partition_geom.dat","r")) == NULL ) {
    fprintf(stderr,"Error:can not open error file!\n");
    return 0;
  }

  for(k=1;k<=kmax;k++)
    for(j=1;j<=jmax;j++)
      for(i=1;i<=imax;i++) {

        fscanf(file_params,"%d" ,&mark); 

        if(mark==1) {
          flagp[FIX(i,j,k)]=1;
          BINDEX[0][index]=i;
          BINDEX[1][index]=j;
          BINDEX[2][index]=k;
          index++;
        }
        delcount++;
        if(delcount==25) {
          fscanf(file_params,"\n"); 
          delcount=0; 
        }
      }
      
  fclose(file_params);
  para->geom->index[1] = index-1;
  return 1;
}

void mark_cell(PARA_DATA *para, REAL **var, int **BINDEX) {
  int i,j, k;
  int ii,ij,ik;
  int SI,EI,SJ,EJ,SK,EK;
  int index=1,pindex=1;
  int imax = para->geom->imax;
  int jmax = para->geom->jmax;
  int kmax = para->geom->kmax;
  int zone_num,zone_inlet,zone_outlet,zone_bl,zone_wall,zone_us,zone_s;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2); 
  REAL *flagp = var[FLAGP],*flagu = var[FLAGU];
  REAL *flagv = var[FLAGV],*flagw = var[FLAGW];
  REAL *flagt = var[FLAGT];

  zone_inlet=para->geom->zone_inlet;
  zone_outlet=para->geom->zone_outlet;
  zone_bl=para->geom->zone_bl;
  zone_wall=para->geom->zone_wall;
  zone_us=para->geom->zone_us;
  zone_s=para->geom->zone_s;

  /**************************************Cell*******************************************************/
  /*
    flag=1 blockage and wall
    flag=0 userdefined
  */
  /**********Blockage**************/
  for(zone_num=zone_outlet+1;zone_num<=zone_bl;zone_num++) {
    SI=BINDEX[7][zone_num];
    EI=BINDEX[8][zone_num];
    SJ=BINDEX[9][zone_num];
    EJ=BINDEX[10][zone_num];
    SK=BINDEX[11][zone_num];
    EK=BINDEX[12][zone_num];
    
    for(ii=SI ;ii<EI ;ii++)
      for(ij=SJ ;ij<EJ ;ij++)
        for(ik=SK ;ik<EK ;ik++) {
          BINDEX[0][index]=ii+1;
          BINDEX[1][index]=ij+1;
          BINDEX[2][index]=ik+1;
          BINDEX[16][FIX(ii+1,ij+1,ik+1)]=zone_num; //zonep
          index++;
          flagp[FIX(ii+1,ij+1,ik+1)]=1; 
        }
  }
  
  
  /**********Premeter Wall Boundary**************/
  ii=0;
  for(ij=0 ;ij<=jmax+1 ;ij++) {
    for(ik=0 ;ik<=kmax+1 ;ik++) {
      BINDEX[0][index]=ii;
      BINDEX[1][index]=ij;
      BINDEX[2][index]=ik;
      BINDEX[16][FIX(ii,ij,ik)]=zone_bl; //zonep
      flagp[FIX(ii,ij,ik)]=1;
      index++;
    }
  }

  ii=imax+1;
  for(ij=0 ;ij<=jmax+1 ;ij++) {
    for(ik=0 ;ik<=kmax+1 ;ik++) {
      BINDEX[0][index]=ii;
      BINDEX[1][index]=ij;
      BINDEX[2][index]=ik;
      BINDEX[16][FIX(ii,ij,ik)]=zone_bl; //zonep
      flagp[FIX(ii,ij,ik)]=1;
      index++;
    }
  }

  ij=0;
  for(ii=0 ;ii<=imax+1 ;ii++) {
    for(ik=0 ;ik<=kmax+1 ;ik++) {
      BINDEX[0][index]=ii;
      BINDEX[1][index]=ij;
      BINDEX[2][index]=ik;
      BINDEX[16][FIX(ii,ij,ik)]=zone_bl; //zonep
      flagp[FIX(ii,ij,ik)]=1;
      index++;
    }
  }

  ij=jmax+1;
  for(ii=0 ;ii<=imax+1 ;ii++) {
    for(ik=0 ;ik<=kmax+1 ;ik++) {
      BINDEX[0][index]=ii;
      BINDEX[1][index]=ij;
      BINDEX[2][index]=ik;
      BINDEX[16][FIX(ii,ij,ik)]=zone_bl; //zonep
      flagp[FIX(ii,ij,ik)]=1;
      index++;
    }
  }

  ik=0;
  for(ii=0 ;ii<=imax+1 ;ii++){
    for(ij=0 ;ij<=jmax+1 ;ij++) {
      BINDEX[0][index]=ii;
      BINDEX[1][index]=ij;
      BINDEX[2][index]=ik;
      BINDEX[16][FIX(ii,ij,ik)]=zone_bl; //zonep
      flagp[FIX(ii,ij,ik)]=1;
      index++;
    }
  }

  ik=kmax+1;
  for(ii=0 ;ii<=imax+1 ;ii++) {
    for(ij=0 ;ij<=jmax+1 ;ij++) {
      BINDEX[0][index]=ii;
      BINDEX[1][index]=ij;
      BINDEX[2][index]=ik;
      BINDEX[16][FIX(ii,ij,ik)]=zone_bl; //zonep
      flagp[FIX(ii,ij,ik)]=1;
      index++;
    }
  }
    
  para->geom->index[1]=index-1;  //block  index for block and premeter wall



/********************************************* User define1d source****************/

  for(zone_num=zone_wall+1;zone_num<=zone_us;zone_num++) {
    SI=BINDEX[7][zone_num];
    EI=BINDEX[8][zone_num];
    SJ=BINDEX[9][zone_num];
    EJ=BINDEX[10][zone_num];
    SK=BINDEX[11][zone_num];
    EK=BINDEX[12][zone_num];
    
    for(ii=SI ;ii<EI ;ii++)
      for(ij=SJ ;ij<EJ ;ij++)
        for(ik=SK ;ik<EK ;ik++) {
          BINDEX[0][index]=ii+1;
          BINDEX[1][index]=ij+1;
          BINDEX[2][index]=ik+1;
          para->geom->iplume=BINDEX[0][index];
          para->geom->jplume=BINDEX[1][index];
          para->geom->kplume=BINDEX[2][index];
          BINDEX[16][FIX(ii+1,ij+1,ik+1)]=zone_num; //zonep
          flagp[FIX(ii+1,ij+1,ik+1)]=0;
          index++;
        }
  }
  
  para->geom->index[2]=index-1;//heat source


  /**************************************Surface boundary**********************/
  /* flag=0 inlet
     flag=1 outlet
     flag=2 wall surface cell
     flag=3 unsigned surface at premeter wall
     flag=4 partition
  */
    
   /**********openings**************/ 
  for(zone_num=zone_us+1;zone_num<=zone_s;zone_num++)  {
    SI=BINDEX[7][zone_num];
    EI=BINDEX[8][zone_num];
    SJ=BINDEX[9][zone_num];
    EJ=BINDEX[10][zone_num];
    SK=BINDEX[11][zone_num];
    EK=BINDEX[12][zone_num];
    if(SI==EI) {
      for(ij=SJ ;ij<EJ ;ij++)
        for(ik=SK ;ik<EK ;ik++) {
          flagu[FIX(SI,ij+1,ik+1)]=-2.0f; // surface cell_openings
        }
    }
    if(SJ==EJ) {
      for(ii=SI ;ii<EI ;ii++)
        for(ik=SK ;ik<EK ;ik++) {
          flagv[FIX(ii+1,SJ,ik+1)]=-2.0f;
        }
    }
    if(SK==EK)  {
      for(ii=SI ;ii<EI ;ii++)
        for(ij=SJ ;ij<EJ ;ij++) {
          flagw[FIX(ii+1,ij+1,SK)]=-2.0f;
        }
    }
  }

    /**********Wall**************/ 
  for(zone_num=zone_bl+1;zone_num<=zone_wall;zone_num++)  {
    SI=BINDEX[7][zone_num];
    EI=BINDEX[8][zone_num];
    SJ=BINDEX[9][zone_num];
    EJ=BINDEX[10][zone_num];
    SK=BINDEX[11][zone_num];
    EK=BINDEX[12][zone_num];
    
    if(SI==EI) {
      for(ij=SJ ;ij<EJ ;ij++)
        for(ik=SK ;ik<EK ;ik++) {
          if(flagu[FIX(SI,ij+1,ik+1)]==-2.0) continue;
          flagu[FIX(SI,ij+1,ik+1)]=2; // surface cell_wall
          BINDEX[13][FIX(SI,ij+1,ik+1)]=zone_num; //zoneu
          
          if(SI<=imax && flagp[FIX(SI,ij+1,ik+1)]<0 && flagp[FIX(SI+1,ij+1,ik+1)]<0){
            flagu[FIX(SI,ij+1,ik+1)]=4; // surface cell-partition
            BINDEX[17][pindex]=SI;
            BINDEX[18][pindex]=ij+1;
            BINDEX[19][pindex]=ik+1;
            pindex++;
          }
        }
    }
  }
  
  para->geom->index[3]=pindex-1;

  for(zone_num=zone_bl+1;zone_num<=zone_wall;zone_num++)  {
    SI=BINDEX[7][zone_num];
    EI=BINDEX[8][zone_num];
    SJ=BINDEX[9][zone_num];
    EJ=BINDEX[10][zone_num];
    SK=BINDEX[11][zone_num];
    EK=BINDEX[12][zone_num];

    if(SJ==EJ) {
      for(ii=SI ;ii<EI ;ii++)
        for(ik=SK ;ik<EK ;ik++) {
          if(flagv[FIX(ii+1,SJ,ik+1)]==-2.0) continue;
          flagv[FIX(ii+1,SJ,ik+1)]=2;
          BINDEX[14][FIX(ii+1,SJ,ik+1)]=zone_num; //zonev
          if(SJ<=jmax && flagp[FIX(ii+1,SJ,ik+1)]<0 && flagp[FIX(ii+1,SJ+1,ik+1)]<0){
            flagv[FIX(ii+1,SJ,ik+1)]=4;
            BINDEX[17][pindex]=ii+1;
            BINDEX[18][pindex]=SJ;
            BINDEX[19][pindex]=ik+1;
            pindex++;
          }
        }
    }
  }
  
  para->geom->index[4]=pindex-1;

  for(zone_num=zone_bl+1;zone_num<=zone_wall;zone_num++)  {
    SI=BINDEX[7][zone_num];
    EI=BINDEX[8][zone_num];
    SJ=BINDEX[9][zone_num];
    EJ=BINDEX[10][zone_num];
    SK=BINDEX[11][zone_num];
    EK=BINDEX[12][zone_num];
    
    if(SK==EK)  {
      for(ii=SI ;ii<EI ;ii++)
        for(ij=SJ ;ij<EJ ;ij++) {
          if(flagw[FIX(ii+1,ij+1,SK)]==-2.0) continue;
          flagw[FIX(ii+1,ij+1,SK)]=2;
          BINDEX[15][FIX(ii+1,ij+1,SK)]=zone_num; //zonev
          if(SK<=kmax && flagp[FIX(ii+1,ij+1,SK)]<0 && flagp[FIX(ii+1,ij+1,SK+1)]<0){
            flagw[FIX(ii+1,ij+1,SK)]=4;
            BINDEX[17][pindex]=ii+1;
            BINDEX[18][pindex]=ij+1;
            BINDEX[19][pindex]=SK;
            pindex++;
          }
        }
    }
  }
 
 
  para->geom->index[5]=pindex-1;

  /****************inlet*******/
  for(zone_num=1;zone_num<=zone_inlet;zone_num++) {
    SI=BINDEX[7][zone_num];
    EI=BINDEX[8][zone_num];
    SJ=BINDEX[9][zone_num];
    EJ=BINDEX[10][zone_num];
    SK=BINDEX[11][zone_num];
    EK=BINDEX[12][zone_num];
    
    if(SI==EI)  {
      for(ij=SJ ;ij<EJ ;ij++)
        for(ik=SK ;ik<EK ;ik++) {
          flagu[FIX(SI,ij+1,ik+1)]=0; //surface cell 0
          BINDEX[13][FIX(SI,ij+1,ik+1)]=zone_num; //zoneu
        }
    }
    
    if(SJ==EJ) {
      for(ii=SI ;ii<EI ;ii++)
        for(ik=SK ;ik<EK ;ik++) {
          flagv[FIX(ii+1,SJ,ik+1)]=0;
          BINDEX[14][FIX(ii+1,SJ,ik+1)]=zone_num; //zonev
        }
    }
    
    if(SK==EK) {
      for(ii=SI ;ii<EI ;ii++)
        for(ij=SJ ;ij<EJ ;ij++) {
          flagw[FIX(ii+1,ij+1,SK)]=0;
          BINDEX[15][FIX(ii+1,ij+1,SK)]=zone_num; //zonew
        }
    }
  }
   
  /****************outlet*******/
  for(zone_num=zone_inlet+1;zone_num<=zone_outlet;zone_num++) {
    SI=BINDEX[7][zone_num];
    EI=BINDEX[8][zone_num];
    SJ=BINDEX[9][zone_num];
    EJ=BINDEX[10][zone_num];
    SK=BINDEX[11][zone_num];
    EK=BINDEX[12][zone_num];

    if(SI==EI) {
      for(ij=SJ ;ij<EJ ;ij++)
        for(ik=SK ;ik<EK ;ik++) {
          flagu[FIX(SI,ij+1,ik+1)]=1; // surface cell_outlet
          BINDEX[13][FIX(SI,ij+1,ik+1)]=zone_num; //zoneu
        }
    }
    
    if(SJ==EJ) {
      for(ii=SI ;ii<EI ;ii++)
        for(ik=SK ;ik<EK ;ik++) {
          flagv[FIX(ii+1,SJ,ik+1)]=1;
          BINDEX[14][FIX(ii+1,SJ,ik+1)]=zone_num; //zoneu
        }
    }
    
    if(SK==EK) {
      for(ii=SI ;ii<EI ;ii++)
        for(ij=SJ ;ij<EJ ;ij++) {
          flagw[FIX(ii+1,ij+1,SK)]=1;
          BINDEX[15][FIX(ii+1,ij+1,SK)]=zone_num; //zoneu
        }
    }
  }
  
  for(index=1;index<=para->geom->index[1];index++) { // All block and premeter 
    i=BINDEX[0][index];
    j=BINDEX[1][index];
    k=BINDEX[2][index];

    if(flagu[FIX(i,j,k)]<0) flagu[FIX(i,j,k)]=3;
    ii=max(i-1,0);  if(flagu[FIX(ii,j,k)]<0) flagu[FIX(ii,j,k)]=3; //block and premeter unsigned face identifier

    if(flagv[FIX(i,j,k)]<0) flagv[FIX(i,j,k)]=3;
    ij=max(j-1,0);  if(flagv[FIX(i,ij,k)]<0) flagv[FIX(i,ij,k)]=3;

    if(flagw[FIX(i,j,k)]<0) flagw[FIX(i,j,k)]=3;
    ik=max(k-1,0);  if(flagw[FIX(i,j,ik)]<0) flagw[FIX(i,j,ik)]=3;
    //if(i==0) printf("%f\t",flagw[FIX(i,j,k)]);
  }

  /****************inlet temperature cells*******/
  for(index=1;index<=para->geom->index[1];index++) { // All solid cells
    i=BINDEX[0][index];
    j=BINDEX[1][index];
    k=BINDEX[2][index];

    if(flagu[FIX(i,j,k)]==0) flagt[FIX(i,j,k)]=1;
    ii=max(i-1,0);  if(flagu[FIX(ii,j,k)]==00) flagt[FIX(i,j,k)]=1; //inlet cells for temperature & other scalars

    if(flagv[FIX(i,j,k)]==0) flagt[FIX(i,j,k)]=1;
    ij=max(j-1,0);  if(flagv[FIX(i,ij,k)]==0) flagt[FIX(i,j,k)]=1;

    if(flagw[FIX(i,j,k)]==0) flagt[FIX(i,j,k)]=1;
    ik=max(k-1,0);  if(flagw[FIX(i,j,ik)]==0) flagt[FIX(i,j,k)]=1;
  }



}





