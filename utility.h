///////////////////////////////////////////////////////////////////////////////
///
/// \file   utility.h
///
/// \brief  Various functions used for data processing in FFD
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
/// This file provides functions for processing the data in FFD
///
///////////////////////////////////////////////////////////////////////////////

REAL check_residual(PARA_DATA *para, REAL **var, REAL *flag,REAL *x);


REAL check_integ(PARA_DATA *para, REAL **var);

int define_IMAX(PARA_DATA *para, int var_type);

int define_JMAX(PARA_DATA *para, int var_type);

void swap(PARA_DATA *para, REAL **var);
void swapT(PARA_DATA *para, REAL **var);
void swapT0(PARA_DATA *para, REAL **var);
void limit(REAL x);

void check_mass(PARA_DATA *para, REAL **var,int i, int j,int k);
REAL check_mass_W(PARA_DATA *para, REAL **var,int i1, int i2,int j1,int j2,int k1,int k2);
REAL check_mass_V(PARA_DATA *para, REAL **var,int i1, int i2,int j1,int j2,int k1,int k2);
REAL check_mass_out(PARA_DATA *para, REAL **var,int i1, int i2,int j1,int j2,int k1,int k2);
REAL check_mass_in(PARA_DATA *para, REAL **var,int i1, int i2,int j1,int j2,int k1,int k2);

void psi_conservation(PARA_DATA *para, REAL **var, REAL *psi,REAL *psi0,int **BINDEX);

REAL outflow(PARA_DATA *para, REAL **var,  REAL *psi,  int **BINDEX);

REAL inflow(PARA_DATA *para, REAL **var,  REAL *psi,  int **BINDEX);

REAL check_max( PARA_DATA *para, REAL *psi, int ci,int cj,int ck);
REAL check_min(PARA_DATA *para, REAL *psi, int ci,int cj,int ck);
REAL check_min_pix(PARA_DATA *para, REAL *psi);
REAL check_max_pix(PARA_DATA *para, REAL *psi);

REAL check_avg(PARA_DATA *para, REAL **var, REAL *psi);
REAL qwall(PARA_DATA *para, REAL **var,int **BINDEX);
REAL check_energy(PARA_DATA *para, REAL **var,int **BINDEX);
void swapuvw(PARA_DATA *para, REAL **var);
void swapuvw_temp(PARA_DATA *para, REAL **var);
void swapUVW_C(PARA_DATA *para, REAL **var);
void swapUVW(PARA_DATA *para, REAL **var);
REAL check_particle(PARA_DATA *para, REAL **var,int i1, int i2,int j1,int j2,int k1,int k2);

REAL num_visco(REAL dx, REAL dt,REAL x_1);
void cfl(PARA_DATA *para, REAL **var);
REAL velmag(REAL u, REAL v);


void u_vortex(PARA_DATA *para, REAL **var, REAL *u, REAL *u0);
void v_vortex(PARA_DATA *para, REAL **var, REAL *v, REAL *v0);
void w_vortex(PARA_DATA *para, REAL **var, REAL *w, REAL *w0);
void p_vortex(PARA_DATA *para, REAL **var, REAL *p, REAL *p0);
void corners(PARA_DATA *para, REAL **var, REAL *psi);
void corners_vortex(PARA_DATA *para, REAL **var, REAL *psi);
void uvw_center(PARA_DATA *para, REAL **var);

REAL qheatsource(PARA_DATA *para, REAL **var,int **BINDEX);
