
void psi_conservation(PARA_DATA *para, REAL **var, int var_type, REAL *psi,REAL *psi0,int **BINDEX);

REAL outflow(PARA_DATA *para, REAL **var,  REAL *psi,  int **BINDEX);

REAL inflow(PARA_DATA *para, REAL **var,  REAL *psi,  int **BINDEX);

REAL check_max( PARA_DATA *para, REAL *psi, int ci,int cj,int ck);
REAL check_min(PARA_DATA *para, REAL *psi, int ci,int cj,int ck);
REAL check_min_pix(PARA_DATA *para, REAL *psi);
REAL check_max_pix(PARA_DATA *para, REAL *psi);

REAL qwall(PARA_DATA *para, REAL **var,int **BINDEX);

void swapuvw(PARA_DATA *para, REAL **var);

void cfl(PARA_DATA *para, REAL **var);

void u_vortex(PARA_DATA *para, REAL **var, REAL *u, REAL *u0);
void v_vortex(PARA_DATA *para, REAL **var, REAL *v, REAL *v0);
void w_vortex(PARA_DATA *para, REAL **var, REAL *w, REAL *w0);
void p_vortex(PARA_DATA *para, REAL **var, REAL *p, REAL *p0);
void corners(PARA_DATA *para, REAL **var, REAL *psi);
void corners_vortex(PARA_DATA *para, REAL **var, REAL *psi);
void uvw_center(PARA_DATA *para, REAL **var);

