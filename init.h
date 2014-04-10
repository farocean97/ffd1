///////////////////////////////////////////////////////////////////////////////
///
/// \file   init.h
///
/// \brief  Subroutines for initializing the FFD
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
/// This file provides functions to initialize the parameters and variables in
/// FFD. The default value for some problem parameters will be assigned.
///
///////////////////////////////////////////////////////////////////////////////

int initial(PARA_DATA *para);

void default_value(PARA_DATA *para);

void free_data (REAL **var);
void free_index(int **BINDEX);

void clear_data(PARA_DATA *para, REAL **var, int **BINDEX);

 