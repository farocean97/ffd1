///////////////////////////////////////////////////////////////////////////////
///
/// \file   solver_gs.h
///
/// \brief  Gauss Seidel equation sovler.
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
/// This file provides functions for solving the discretize equation by
/// iterations. In FFD, the diffusion equation and pressure equaion are solved
/// through iterative approaches.
///
///////////////////////////////////////////////////////////////////////////////

void Gauss_Seidel(PARA_DATA *para, REAL **var, REAL *flag, REAL *x);

void GS_P(PARA_DATA *para, REAL **var, int Type, REAL *x);
void Gauss_Seidel_simple(PARA_DATA *para, REAL **var, REAL *flag, REAL *x);