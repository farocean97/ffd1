///////////////////////////////////////////////////////////////////////////////
///
/// \file   init.c
///
/// \brief  Subroutines for initializing the FFD
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
/// This file provides functions to initialize the parameters and variables in
/// FFD. The default value for some problem parameters will be assigned.
///
///////////////////////////////////////////////////////////////////////////////
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "data_structure.h"
#include "init.h"
#include "parameters.h"


///////////////////////////////////////////////////////////////////////////////
///\brief Entrance of initializing the parameter
///
///\param para Pointer to FFD parameters
///
///\return 1 if no error occurs
///////////////////////////////////////////////////////////////////////////////
int initial(PARA_DATA *para) {

  default_value(para);
  input_para(para);

  return(1);
} // End of initial( )

///////////////////////////////////////////////////////////////////////////////
///\brief Specifying default value for problem parameters.
///
///\param para Pointer to FFD parameters
///
///\return void No return needed
///////////////////////////////////////////////////////////////////////////////
void default_value(PARA_DATA *para){
  para->mytime->dt = 0.1f; 
  para->mytime->t_output = 1000;
  para->mytime->t  = 0.0;
  para->mytime->t_start = 0.0;
  para->mytime->t_step = 0;
  para->mytime->t_start = clock();


  para->prob->chen_a = 0.03874;  /* the coeffcient of Chen's model*/
  para->prob->Prt = 0.9;   /* turbulent Prandl number */ 
   para->prob->mu = 0.01;  /* mu_air */
  para->prob->rho = 1.0;
  para->prob->tur_model = LAM;
  para->bc->qs=0;

  para->solv->check_residual = 0; /* donot check residual */
  para->solv->solver = GS;        /* Gauss-Seidel Solver */

  para->outp->Temp_ref   = 300.0;//35.5f;//10.25f;

} // End of default_value()

///////////////////////////////////////////////////////////////////////////////
///\brief Freeing the memeory allocated by the variables
///
///\param var Pointer to FFD variables
///
///\return void No return needed
///////////////////////////////////////////////////////////////////////////////
void free_data(REAL **var) {
  if(var[X]) free(var[X]);
  if(var[Y]) free(var[Y]);
  if(var[Z]) free(var[Z]);
  if(var[VX]) free(var[VX]);
  if(var[VY]) free(var[VY]);
  if(var[VZ]) free(var[VZ]);
  if(var[VXS]) free(var[VXS]);
  if(var[VYS]) free(var[VYS]);
  if(var[VZS]) free(var[VZS]);
  if(var[VXM]) free(var[VXM]);
  if(var[VYM]) free(var[VYM]);
  if(var[VZM]) free(var[VZM]);
  if(var[DEN]) free(var[DEN]);
  if(var[DENS]) free(var[DENS]);
  if(var[TEMP]) free(var[TEMP]);
  if(var[TEMPM]) free(var[TEMPM]);
  if(var[TEMPS]) free(var[TEMPS]);
  if(var[IP]) free(var[IP]);
  if(var[TMP1]) free(var[TMP1]);
  if(var[TMP2]) free(var[TMP2]);
  if(var[TMP3]) free(var[TMP3]);
  if(var[AP]) free(var[AP]);
  if(var[AN]) free(var[AN]);
  if(var[AS]) free(var[AS]);
  if(var[AE]) free(var[AE]);
  if(var[AW]) free(var[AW]);
  if(var[AF]) free(var[AF]);
  if(var[AB]) free(var[AB]);
  if(var[B])  free(var[B]);
  if(var[GX])  free(var[GX]);
  if(var[GY])  free(var[GY]);
  if(var[GZ])  free(var[GZ]);
  if(var[AP0])  free(var[AP0]);
  if(var[PP])  free(var[PP]);
  if(var[FLAGP])  free(var[FLAGP]);
  if(var[FLAGU])  free(var[FLAGU]);
  if(var[FLAGV])  free(var[FLAGV]);
  if(var[FLAGW])  free(var[FLAGW]);
  if(var[LOCMIN])  free(var[LOCMIN]);
  if(var[LOCMAX])  free(var[LOCMAX]);
  if(var[VXBC])  free(var[VXBC]);
  if(var[VYBC])  free(var[VYBC]);
  if(var[VZBC])  free(var[VZBC]);
  if(var[TEMPBC])  free(var[TEMPBC]);
  if(var[FX])  free(var[FX]);
  if(var[FY])  free(var[FY]);
  if(var[FZ])  free(var[FZ]);
  if(var[ER])  free(var[ER]);
  if(var[RR])  free(var[RR]);
  if(var[VXC]) free(var[VXC]);
  if(var[VYC]) free(var[VYC]);
  if(var[VZC]) free(var[VZC]);
   if(var[FLUX]) free(var[FLUX]);
  if(var[VSX]) free(var[VSX]);
  if(var[VSY]) free(var[VSY]);
  if(var[VSZ]) free(var[VSZ]);
  if(var[CFLX]) free(var[CFLX]);
  if(var[CFLY]) free(var[CFLY]);
  if(var[CFLZ]) free(var[CFLZ]);
  if(var[DIST]) free(var[DIST]);
  if(var[VT]) free(var[VT]);
  if(var[FLAGT]) free(var[FLAGT]);
} // End free_data()


///////////////////////////////////////////////////////////////////////////////
///\brief Freeing the memeory of boundary indices
///
///\param BINDEX Pointer to FFD bounary indices
///
///\return void No return needed
///////////////////////////////////////////////////////////////////////////////
void free_index(int **BINDEX) {
  if(BINDEX[0]) free(BINDEX[0]);
  if(BINDEX[1]) free(BINDEX[1]);
  if(BINDEX[2]) free(BINDEX[2]);
  if(BINDEX[3]) free(BINDEX[3]);
  if(BINDEX[4]) free(BINDEX[4]);
  if(BINDEX[5]) free(BINDEX[5]);
  if(BINDEX[6]) free(BINDEX[6]);
  if(BINDEX[7]) free(BINDEX[7]);
  if(BINDEX[8]) free(BINDEX[8]);
  if(BINDEX[9]) free(BINDEX[9]);
  if(BINDEX[10]) free(BINDEX[10]);
  if(BINDEX[11]) free(BINDEX[11]);
  if(BINDEX[12]) free(BINDEX[12]);
  if(BINDEX[13]) free(BINDEX[13]);
  if(BINDEX[14]) free(BINDEX[14]);
  if(BINDEX[15]) free(BINDEX[15]);
  if(BINDEX[16]) free(BINDEX[16]);
  if(BINDEX[17]) free(BINDEX[17]);
  if(BINDEX[18]) free(BINDEX[18]);
  if(BINDEX[19]) free(BINDEX[19]);
  if(BINDEX[20]) free(BINDEX[20]);
}

///////////////////////////////////////////////////////////////////////////////
///\brief Initializing value for FFD variables
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD variables
///\param BINDEX Pointer to FFD boundary index
///
///\return void No return needed
///////////////////////////////////////////////////////////////////////////////
void clear_data(PARA_DATA *para, REAL **var, int **BINDEX) {
  int i, size;
  size=(para->geom->imax+2) * (para->geom->jmax+2) * (para->geom->kmax+2);
  
  para->mytime->t = 0.0;
  para->mytime->t_step = 0;

  for(i=0; i<size; i++)  {
    var[GX][i]     = 0.0;
    var[GY][i]     = 0.0;
    var[GZ][i]     = 0.0;
    var[VX][i]     = 0.0;
    var[VY][i]     = 0.0;
    var[VZ][i]     = 0.0;
    var[VXM][i]    = 0.0;
    var[VYM][i]    = 0.0;
    var[VZM][i]    = 0.0;
    var[VXC][i]    = 0.0;
    var[VYC][i]    = 0.0;
    var[VZC][i]    = 0.0;
    var[VXS][i]    = 0.0;
    var[VYS][i]    = 0.0;
    var[VZS][i]    = 0.0;
    var[DEN][i]    = 0.0;    
    var[DENS][i]   = 0.0;
    var[TEMP][i]   = 293.0;
    var[TEMPM][i]  = 0.0;
    var[TEMPS][i]  = 0.0;
    var[IP][i]     = 0.0;
    var[AP][i]     = 0.0;
    var[AW][i]     = 0.0;
    var[AE][i]     = 0.0;
    var[AS][i]     = 0.0;
    var[AN][i]     = 0.0;
    var[AB][i]     = 0.0;
    var[AF][i]     = 0.0;
    var[B][i]      = 0.0;
    var[AP0][i]    = 0.0;
    var[TMP1][i]   = 0.0;
    var[TMP2][i]   = 0.0;
    var[TMP3][i]   = 0.0;
    var[PP][i]     = 0.0;
    var[FLAGP][i]    = -1.0;
    var[FLAGU][i]    = -1.0;
    var[FLAGV][i]    = -1.0;
    var[FLAGW][i]    = -1.0;
    var[FLAGT][i]    = -1.0;
    var[VXBC][i]     = 0.0;
    var[VYBC][i]     = 0.0;
    var[VZBC][i]     = 0.0;
    var[TEMPBC][i]   = 0.0;
    var[ER][i]       =0.0;
    var[RR][i]       =0.0;
    var[FLUX][i]     =0.0;
    var[VSX][i]     = 0.0;
    var[VSY][i]     = 0.0;
    var[VSZ][i]     = 0.0;
    var[CFLX][i]     = 0.0;
    var[CFLY][i]     = 0.0;
    var[CFLZ][i]     = 0.0;
    var[DIST][i]     =0.0;
    var[VT][i]       =0.0;
    BINDEX[20][i]    =0.0;
  }

} // End of clear_data()

