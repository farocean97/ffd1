/* ----------------------------------------------------------------------------

  File:          write.h

  Written by:   Wangda Zuo

  Task:         Declare subroutines in write.c              

---------------------------------------------------------------------------- */

int write_data(PARA_DATA *para, REAL **var, char *name,int **BINDEX);

int write_data1(PARA_DATA *para, REAL **var, char *name);

int write_unsteady(PARA_DATA *para, REAL **var, char *name);

int write_SCI(PARA_DATA *para, REAL **var, char *name);

int write_time(int tsize, REAL *T_mon, REAL *U_mon,REAL *P_mon);

int write_displacement(PARA_DATA *para, REAL **var);
int write_data_block(PARA_DATA *para, REAL **var, char *name);
void write_temperature(PARA_DATA *para, REAL **var, int **BINDEX);

void write_plume_mass(PARA_DATA *para, REAL **var, int **BINDEX);
