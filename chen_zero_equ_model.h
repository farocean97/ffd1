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

REAL nu_t_chen_zero_equ(REAL u,REAL v,REAL w, REAL dist);
void lengthscale(PARA_DATA *para, REAL **var);
void vis_turb(PARA_DATA *para, REAL **var);
