void advection(PARA_DATA *para, REAL **var, int var_type,
                REAL *flag, REAL *d, REAL *d0, int **BINDEX);

void traceback(PARA_DATA *para, REAL **var, int **BINDEX);
void traceback_UVW(PARA_DATA *para, REAL **var, int var_type, REAL *flag, int **BINDEX);



void XLOCATION(PARA_DATA *para, REAL **var,  REAL *flag, 
               REAL *x, REAL u0, int i, int j, int k,  
               REAL *OL, int *OC, int *LOC , int *COOD);
void YLOCATION(PARA_DATA *para, REAL **var,  REAL *flag,
               REAL *y, REAL v0, int i, int j, int k,  
               REAL *OL, int *OC, int *LOC , int *COOD);
void ZLOCATION(PARA_DATA *para, REAL **var,  REAL *flag, 
               REAL *z, REAL w0, int i, int j, int k,  
               REAL *OL, int *OC, int *LOC , int *COOD);
void XLOCATION_U(PARA_DATA *para, REAL **var, REAL *flag, REAL *x, REAL u0,
               int i, int j, int k,  REAL *OL, int *OC, int *LOC , int *COOD);
void YLOCATION_V(PARA_DATA *para, REAL **var, REAL *flag, REAL *y, REAL v0,
               int i, int j, int k,  REAL *OL, int *OC, int *LOC , int *COOD);
void ZLOCATION_W(PARA_DATA *para, REAL **var, REAL *flag, REAL *z, REAL w0,
               int i, int j, int k,  REAL *OL, int *OC, int *LOC , int *COOD);