///////////////////////////////////////////////////////////////////////////////
///
/// \file   chen_zero_equ_model.h
///
/// \brief  Subroutines for zeo equation turbulence model
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
/// This file provides functions that are used for the computing turbulence
/// viscosity based on the zero equation model.
///////////////////////////////////////////////////////////////////////////////
REAL nu_t_chen_zero_equ(REAL u,REAL v,REAL w, REAL dist);
void lengthscale(PARA_DATA *para, REAL **var);
void vis_turb(PARA_DATA *para, REAL **var);
