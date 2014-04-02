///////////////////////////////////////////////////////////////////////////////
///
/// \file   advection.h
///
/// \brief  Solver for advection step
///
/// \author Mingang Jin, Qingyan Chen
///         Purdue University
///         Jin56@purdue.edu, YanChen@purdue.edu
///         Wangda Zuo
///         University of Miami
///         W.Zuo@miami.edu
///
/// \date   8/3/2013
///
/// This file provides functions that used for the advection step of FFD method.
/// The advection starts with \c advect(). Then different subroutines are 
/// called according to the properties of the variables that are sorted by
/// the location of variables assigned in the control volume. 
/// Velocities at X, Y and Z directions are locatted
/// on the surface of the control volume. They are computed using 
/// subroutines: \c traceback_uvw
/// Scalar variables are in the center of control volume and they are computed
/// using \c trace_scalar().
///
///////////////////////////////////////////////////////////////////////////////

void traceback(PARA_DATA *para, REAL **var, int **BINDEX);
void traceback_UVW(PARA_DATA *para, REAL **var, int var_type, REAL *flag, int **BINDEX);

void advection(PARA_DATA *para, REAL **var, int var_type,
                REAL *flag, REAL *d, REAL *d0, int **BINDEX);

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