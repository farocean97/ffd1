///////////////////////////////////////////////////////////////////////////////
///
/// \file   write_data.h
///
/// \brief  Subroutines for exporting the simulation results
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
/// This file provides functions for exporting the simulation results.FFD exports
/// results that can be read by SCI through \c write_SCI(),and also exports data
/// format that can be open by Tecplot through \c write_Tecplot(). In addition,
/// data without any treatmen was exported by \c write_unsteady. This data file 
/// will be used, if FFD continues the simulation. 
///
///////////////////////////////////////////////////////////////////////////////

int write_Tecplot(PARA_DATA *para, REAL **var, char *name,int **BINDEX);

int write_unsteady(PARA_DATA *para, REAL **var, char *name);

int write_SCI(PARA_DATA *para, REAL **var, char *name);

