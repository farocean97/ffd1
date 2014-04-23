///////////////////////////////////////////////////////////////////////////////
///
/// \file   read_data.c
///
/// \brief  Reading data from previous simulation results
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
/// This file provides functions for reading date from previous simualtion results
/// so that FFD can continue the simulation. 
///
///////////////////////////////////////////////////////////////////////////////
#define _CRT_SECURE_NO_DEPRECATE
#include <stdio.h>
#include <string.h>

#include "data_structure.h"
#include "read_data.h"

FILE *file_params;

///////////////////////////////////////////////////////////////////////////////
///\brief Read the unsteady result and continues simulation
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///
///\return 1 no error
///////////////////////////////////////////////////////////////////////////////
int read_data(PARA_DATA *para, REAL **var) {
  int i,j, k;
  int imax = para->geom->imax;
  int jmax = para->geom->jmax;
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  char string[400];


  if( (file_params=fopen("unsteady.plt","r")) == NULL ) {
    fprintf(stderr,"Error:can not open error file!\n");
    return 0;
  }
 
  FOR_ALL_CELL
    fgets(string, 400, file_params);
    sscanf( string,"%f%f%f%f%f%f", &var[VX][FIX(i,j,k)], &var[VY][FIX(i,j,k)],
                                   &var[VZ][FIX(i,j,k)], &var[TEMP][FIX(i,j,k)],
                                   &var[DEN][FIX(i,j,k)], &var[IP][FIX(i,j,k)]);
  END_FOR
    
  fclose(file_params);

  return 1;

} // End of read_dara()

