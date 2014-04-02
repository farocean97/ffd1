#define _CRT_SECURE_NO_DEPRECATE
#include <stdio.h>
#include <string.h>

#include "data_structure.h"
#include "read_data.h"

FILE *file_params;

/******************************************************************************
| Write the data to a file for Tecplot 
******************************************************************************/
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
    sscanf( string,"%f%f%f%f%f%f", &var[VX][FIX(i,j,k)], &var[VY][FIX(i,j,k)], &var[VZ][FIX(i,j,k)], &var[TEMP][FIX(i,j,k)],
                                   &var[DEN][FIX(i,j,k)], &var[IP][FIX(i,j,k)]);
  END_FOR
    
  fclose(file_params);

  return 1;

} // End of read_dara()

