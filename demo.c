/*======================================================================
 demo.c
=======================================================================*/
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#include "data_structure.h"
#include "solver.h"
#include "write_data.h"
#include "init.h"
#include "boundary.h"
#include "read_data.h"

#include "timing.h"
#include "input.h"

/* global variables */
static REAL dt, diff, visc;
static REAL force, source;
static int screen;

REAL **var;
int  **BINDEX;
int  *xindex,*yindex,*zindex,*fltemp;
int  *coodx,*coody,*coodz;
int  *zonex1,*zonex2,*zoney1,*zoney2,*zonez1,*zonez2;
int  *zoneu,*zonev,*zonew,*zonep,*pzonex,*pzoney,*pzonez;
int  *plumdex;
REAL *x, *y, *z, *gx, *gy, *gz;
REAL *u, *v, *w, *u_s, *v_s, *w_s, *u_mean, *v_mean, *w_mean;
REAL *uc,*vc,*wc;
REAL *dens, *dens_s, *temp, *temp_s, *temp_mean, *p, *my_div, *pp;
REAL *tmp1, *tmp2, *tmp3;
REAL *ap, *an, *as, *aw, *ae, *b, *ab, *af, *ap0;
REAL *flagp, *flagu, *flagv, *flagw,*flagt;
REAL *locmin,*locmax;
REAL *vxbc,*vybc,*vzbc,*tempbc;
REAL *fx,*fy,*fz;
REAL *er,*rr;
REAL *flux;
REAL *vsx,*vsy,*vsz;
REAL *cflx,*cfly,*cflz;
REAL *dist,*vt;

static GEOM_DATA geom;
static PROB_DATA prob;
static TIME_DATA mytime;
static OUTP_DATA outp1;
static BC_DATA bc;
static SOLV_DATA solv;
static PARA_DATA para;


clock_t start, end;

/******************************************************************************
| allocate data
******************************************************************************/
int allocate_data ( void ){
  int size = (geom.imax+2) * (geom.jmax+2) * (geom.kmax+2);
  printf( "size=%d\n", size); 
  x          = (REAL *) malloc ( size*sizeof(REAL) );
  y          = (REAL *) malloc ( size*sizeof(REAL) );
  z          = (REAL *) malloc ( size*sizeof(REAL) );
  u          = (REAL *) malloc ( size*sizeof(REAL) );
  v          = (REAL *) malloc ( size*sizeof(REAL) );
  w          = (REAL *) malloc ( size*sizeof(REAL) );
  u_s        = (REAL *) malloc ( size*sizeof(REAL) );
  v_s        = (REAL *) malloc ( size*sizeof(REAL) );
  w_s        = (REAL *) malloc ( size*sizeof(REAL) );
  u_mean    = (REAL *) malloc ( size*sizeof(REAL) );
  v_mean    = (REAL *) malloc ( size*sizeof(REAL) );
  w_mean    = (REAL *) malloc ( size*sizeof(REAL) );
  temp      = (REAL *) malloc ( size*sizeof(REAL) );
  temp_s     = (REAL *) malloc ( size*sizeof(REAL) );
  temp_mean  = (REAL *) malloc ( size*sizeof(REAL) );
  dens      = (REAL *) malloc ( size*sizeof(REAL) );  
  dens_s    = (REAL *) malloc ( size*sizeof(REAL) );
  p         = (REAL *) malloc ( size*sizeof(REAL) ); 
  tmp1      = (REAL *) malloc ( size*sizeof(REAL) );  
  tmp2      = (REAL *) malloc ( size*sizeof(REAL) );  
  tmp3      = (REAL *) malloc ( size*sizeof(REAL) );  
  ap        = (REAL *) malloc ( size*sizeof(REAL) );
  an        = (REAL *) malloc ( size*sizeof(REAL) );
  as        = (REAL *) malloc ( size*sizeof(REAL) );
  aw        = (REAL *) malloc ( size*sizeof(REAL) );
  ae        = (REAL *) malloc ( size*sizeof(REAL) );
  ab        = (REAL *) malloc ( size*sizeof(REAL) );
  af        = (REAL *) malloc ( size*sizeof(REAL) );
  b        = (REAL *) malloc ( size*sizeof(REAL) );
  gx         = (REAL *) malloc ( size*sizeof(REAL) );  
  gy       = (REAL *) malloc ( size*sizeof(REAL) );
  gz       = (REAL *) malloc ( size*sizeof(REAL) );
  ap0      = (REAL *) malloc ( size*sizeof(REAL) );
  pp      = (REAL *) malloc ( size*sizeof(REAL) );
  flagp      = (REAL *) malloc ( size*sizeof(REAL) );
  flagu      = (REAL *) malloc ( size*sizeof(REAL) );
  flagv      = (REAL *) malloc ( size*sizeof(REAL) );
  flagw      = (REAL *) malloc ( size*sizeof(REAL) );
  locmin      = (REAL *) malloc ( size*sizeof(REAL) );
  locmax      = (REAL *) malloc ( size*sizeof(REAL) );
  vxbc      = (REAL *) malloc ( size*sizeof(REAL) );
  vybc      = (REAL *) malloc ( size*sizeof(REAL) );
  vzbc      = (REAL *) malloc ( size*sizeof(REAL) );
  tempbc      = (REAL *) malloc ( size*sizeof(REAL) );
  fx      = (REAL *) malloc ( size*sizeof(REAL) );
  fy      = (REAL *) malloc ( size*sizeof(REAL) );
  fz      = (REAL *) malloc ( size*sizeof(REAL) );
  er      = (REAL *) malloc ( size*sizeof(REAL) ); 
  rr      = (REAL *) malloc ( size*sizeof(REAL) ); 
  uc          = (REAL *) malloc ( size*sizeof(REAL) );
  vc          = (REAL *) malloc ( size*sizeof(REAL) );
  wc          = (REAL *) malloc ( size*sizeof(REAL) );
  flux       = (REAL *) malloc ( size*sizeof(REAL) );
  vsx          = (REAL *) malloc ( size*sizeof(REAL) );
  vsy          = (REAL *) malloc ( size*sizeof(REAL) );
  vsz          = (REAL *) malloc ( size*sizeof(REAL) );
  cflx          = (REAL *) malloc ( size*sizeof(REAL) );
  cfly          = (REAL *) malloc ( size*sizeof(REAL) );
  cflz          = (REAL *) malloc ( size*sizeof(REAL) );
  dist          = (REAL *) malloc ( size*sizeof(REAL) );
  vt          = (REAL *) malloc ( size*sizeof(REAL) );
  flagt          = (REAL *) malloc ( size*sizeof(REAL) );
  var        = (REAL **) malloc ( 62*sizeof(REAL*) );
  
  var[X]     = x;
  var[Y]     = y;
  var[Z]     = z;
  var[VX]    = u;
  var[VY]    = v;
  var[VZ]    = w;
  var[VXS]   = u_s;
  var[VYS]   = v_s;
  var[VZS]   = w_s;
  var[VXM]   = u_mean;
  var[VYM]   = v_mean;
  var[VZM]   = w_mean;
  var[DEN]   = dens;
  var[DENS]  = dens_s;
  var[IP]    = p;
  var[TEMP]  = temp;
  var[TEMPS] = temp_s;
  var[TEMPM] = temp_mean;
  var[AP]    = ap;
  var[AN]    = an;
  var[AS]    = as;
  var[AW]    = aw;
  var[AE]    = ae;
  var[AB]    = ab;
  var[AF]    = af;
  var[B]     = b;
  var[TMP1]  = tmp1;
  var[TMP2]  = tmp2;
  var[TMP3]  = tmp3;
  var[GX]    = gx;
  var[GY]    = gy;
  var[GZ]    = gz;
  var[AP0]   = ap0;
  var[PP]    = pp;
  var[FLAGP] =flagp;
  var[FLAGU] =flagu;
  var[FLAGV] =flagv;
  var[FLAGW] =flagw;
  var[LOCMIN] =locmin;
  var[LOCMAX] =locmax;
  var[VXBC] =vxbc;
  var[VYBC] =vybc;
  var[VZBC] =vzbc;
  var[TEMPBC] =tempbc;
  var[FX] =fx;
  var[FY] =fy;
  var[FZ] =fz;
  var[ER] =er;
  var[RR] =rr;
  var[VXC]    = uc;
  var[VYC]    = vc;
  var[VZC]    = wc;
  var[FLUX]   =flux;
  var[VSX]    = vsx;
  var[VSY]    = vsy;
  var[VSZ]    = vsz;
  var[CFLX]    = cflx;
  var[CFLY]    = cfly;
  var[CFLZ]    = cflz;
  var[DIST] =dist;
  var[VT]   =vt;
  var[FLAGT]   =flagt;

  xindex          = (int *) malloc ( size*sizeof(int) );
  yindex          = (int *) malloc ( size*sizeof(int) );
  zindex          = (int *) malloc ( size*sizeof(int) );
  fltemp          = (int *) malloc ( size*sizeof(int) );
  coodx              = (int *) malloc ( size*sizeof(int) );
  coody            = (int *) malloc ( size*sizeof(int) );
  coodz           = (int *) malloc ( size*sizeof(int) );
  zonex1           = (int *) malloc ( size*sizeof(int) );
  zonex2           = (int *) malloc ( size*sizeof(int) );
  zoney1           = (int *) malloc ( size*sizeof(int) );
  zoney2           = (int *) malloc ( size*sizeof(int) );
  zonez1           = (int *) malloc ( size*sizeof(int) );
  zonez2           = (int *) malloc ( size*sizeof(int) );
  zoneu           = (int *) malloc ( size*sizeof(int) );
  zonev           = (int *) malloc ( size*sizeof(int) );
  zonew           = (int *) malloc ( size*sizeof(int) );
  zonep           = (int *) malloc ( size*sizeof(int) );
  pzonex           = (int *) malloc ( size*sizeof(int) );
  pzoney           = (int *) malloc ( size*sizeof(int) );
  pzonez           = (int *) malloc ( size*sizeof(int) );
  plumdex          = (int *) malloc ( size*sizeof(int) );

  BINDEX              = (int **) malloc ( 21*sizeof(int*) );

  BINDEX[0]= xindex;
  BINDEX[1]= yindex;
  BINDEX[2]= zindex;
  BINDEX[3]= fltemp;
  BINDEX[4]= coodx;
  BINDEX[5]= coody;
  BINDEX[6]= coodz;
  BINDEX[7]= zonex1;
  BINDEX[8]= zonex2;
  BINDEX[9]= zoney1;
  BINDEX[10]= zoney2;
  BINDEX[11]= zonez1;
  BINDEX[12]= zonez2;
  BINDEX[13]= zoneu;
  BINDEX[14]= zonev;
  BINDEX[15]= zonew;
  BINDEX[16]= zonep;
  BINDEX[17]= pzonex;
  BINDEX[18]= pzoney;
  BINDEX[19]= pzonez;
  BINDEX[20]= plumdex;

  if( !x || !y || !z || !u || !v || !w || !u_s || !v_s || !w_s || 
      !u_mean || !v_mean || !w_mean || 
      !dens || !dens_s || !temp || !temp_s || !temp_mean || 
      !tmp1 || !tmp2 || !tmp3 ||
      !ap || !ae || !aw || !as || !an || !ab || !af || !b || !gx || !gy || !gz || !ap0 || !pp || !flagp ||
      ! flagu || ! flagv || ! flagw || ! locmin || ! locmax ||
      ! vxbc ||! vybc ||! vzbc ||! tempbc||! xindex ||! yindex ||! zindex||! fx||! fy||! fz|| !er || !rr ||
      !uc || !vc || !wc || !flux|| !vsx|| !vsy|| !vsz|| !cflx|| !cfly|| !cflz ||!dist || !vt|| !flagt) {
    fprintf ( stderr, "cannot allocate data\n" );
    return ( 0 );
  }
  return ( 1 );
} /** allocate_data() **/


/******************************************************************************
   main --- main routine
******************************************************************************/
int main() { 
  para.geom = &geom;
  para.outp = &outp1;
  para.prob = &prob;
  para.mytime = &mytime;
  para.bc     = &bc;
  para.solv   = &solv;

  printf("hello");
  initial(&para);

  printf("hello");
  
  if(!read_max(&para, var)) {
    printf("no file"); exit(1);
  }

  //printf("imax= %d\t jmax= %d\t  kmax= %d\n ", para.geom->imax,para.geom->jmax,para.geom->kmax);

  if(!allocate_data( ))    exit ( 1 );

  clear_data(&para,  var,BINDEX);
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
  if(!read_input(&para, var,BINDEX)) {
    printf("no file"); exit(1);
  }

 // if(!read_zeroone(&para, var,BINDEX)) {printf("no file"); exit(1);}

  mark_cell(&para, var,BINDEX);

  if(para.solv->read_file==1) read_data(&para,var);
  
  FFD_solver(&para, var,BINDEX);
  write_SCI(&para,  var, "output");
  free_data(var);
  free_index(BINDEX);
  
  getchar();  
  exit ( 0 );
} // End of main( )