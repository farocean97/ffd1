///////////////////////////////////////////////////////////////////////////////
///
/// \file   bounary.h
///
/// \brief  Subroutines for setting boundary conditions
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
/// This file provides functions that used for the setting boundary for FFD.
/// FFD will set up bounary conditons for velocity \c set_vel_bnd(), pressure
/// \c set_bnd_pressure(),temperature \c set_bnd_temp() and species 
/// concentration \c set_bnd_density(). In addition, this fine contains 
/// functions for enforcing mass conservation at the outflow bounary and 
/// special treatment for plume cells.
///
///////////////////////////////////////////////////////////////////////////////

//void set_vel_bnd(PARA_DATA *para, REAL **var,int **BINDEX);

void set_bnd(PARA_DATA *str_geom, REAL **var, int var_type, REAL *x, int **BINDEX);

void set_bnd_temp(PARA_DATA *para, REAL **var, int var_type, REAL *psi,int **BINDEX);

void set_bnd_pressure(PARA_DATA *para, REAL **var, REAL *p, int **BINDEX);

void set_bnd_density(PARA_DATA *para, REAL **var, int var_type, REAL *p,int **BINDEX);

void set_bnd_vel(PARA_DATA *para, REAL **var, int var_type, REAL *vx, int **BINDEX);

void mass_conservation(PARA_DATA *para, REAL **var,int **BINDEX);

void plume(PARA_DATA *para, REAL **var, int **BINDEX);

void plume_thermal(PARA_DATA *para, REAL **var, int **BINDEX);
