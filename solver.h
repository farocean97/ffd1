///////////////////////////////////////////////////////////////////////////////
///
/// \file   solver.h
///
/// \brief  The solver of FFD
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
/// This file provides functions for performing the computation of solving the
/// governing equation.FFD solves the set of equations sequentially.First the
/// \c vel_step() is applied to solve the momentum equation, and then energy 
/// equation \c temp_step(), and species transport \c den_step().
///
///////////////////////////////////////////////////////////////////////////////



void den_step(PARA_DATA *para, REAL **var,int **BINDEX);

void vel_step(PARA_DATA *para, REAL **var, int **BINDEX);

void temp_step(PARA_DATA *para, REAL **var,int **BINDEX);

void FFD_solver(PARA_DATA *para, REAL **var, int **BINDEX);

void equ_solver(PARA_DATA *para, REAL **var, int Type, REAL *x);
