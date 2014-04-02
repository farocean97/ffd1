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
/// This file provides functions that used for the setting boundary for FFD.
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

void set_vel_bnd(PARA_DATA *para, REAL **var,int **BINDEX);

void set_bnd(PARA_DATA *str_geom, REAL **var, int var_type, REAL *x, int **BINDEX);

void set_bnd_temp(PARA_DATA *para, REAL **var, int var_type, REAL *psi,int **BINDEX);

void set_bnd_pressure(PARA_DATA *para, REAL **var, REAL *p, int **BINDEX);

void set_bnd_density(PARA_DATA *para, REAL **var, int var_type, REAL *p,int **BINDEX);

void set_bnd_vel(PARA_DATA *para, REAL **var, int var_type, REAL *vx, int **BINDEX);

void mass_conservation(PARA_DATA *para, REAL **var,int **BINDEX);

void plume(PARA_DATA *para, REAL **var, int **BINDEX);

void plume_thermal(PARA_DATA *para, REAL **var, int **BINDEX);
