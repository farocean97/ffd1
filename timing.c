///////////////////////////////////////////////////////////////////////////////
///
/// \file   timing.c
///
/// \brief  Counting the computing time.
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
/// This file provides functions for counting the computing time. The CPU time
/// and the physical time are recorded during the simulaiton for comparison.
///
///////////////////////////////////////////////////////////////////////////////

#include <stdio.h> 
#include <time.h>

#include "data_structure.h"
#include "timing.h"

///////////////////////////////////////////////////////////////////////////////
///\brief Calculate the simulation time and time ratio
///
///\param para Pointer to FFD parameters
///
///\return No return needed
///////////////////////////////////////////////////////////////////////////////
void timing(PARA_DATA *para) {
  double cputime;

  para->mytime->t += para->mytime->dt;
  para->mytime->t_step += 1;
  para->mytime->t_end = clock();

  cputime= ((double) (clock() - para->mytime->t_start) / CLOCKS_PER_SEC);

  printf("Phyical time=%.1f s, CPU time=%.3f s, Time Ratio=%.4f\n", 
  para->mytime->t, cputime, para->mytime->t/cputime);

} // End of timing( )