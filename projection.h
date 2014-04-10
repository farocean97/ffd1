///////////////////////////////////////////////////////////////////////////////
///
/// \file   projection.h
///
/// \brief  Subroutines for performing pressure projection used in FFD
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
/// This file provides functions for performing pressure projection in FFD, which
/// is realized in two steps,including solving the pressure equation and performing
/// the velocity correction.
///
///////////////////////////////////////////////////////////////////////////////

void projection(PARA_DATA *para, REAL **var, int **BINDEX);
void vel_correction(PARA_DATA *para, REAL **var,int **BINDEX);