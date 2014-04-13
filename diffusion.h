///////////////////////////////////////////////////////////////////////////////
///
/// \file   diffusion.h
///
/// \brief  Subroutines for solving the diffusion equation
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
/// This file provides functions for solving the diffusion equation \c diffusion().
/// The coefficients for the discretized equation will be first assigned and then
/// it will be sovled by iterative approach.
///
///////////////////////////////////////////////////////////////////////////////

void diffusion(PARA_DATA *para, REAL **var, int var_type, REAL *psi, REAL *psi0,int **BINDEX);

void coef_diff(PARA_DATA *para, REAL **var, REAL *psi, REAL *psi0, int var_type, int **BINDEX);

void source(PARA_DATA *para, REAL **var, int var_type,  REAL *psi, REAL *psi0, int **BINDEX);