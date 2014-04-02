#include "data_structure.h"

/******************************************************************************
|  Input Parameters
******************************************************************************/
void input_para(PARA_DATA *para) {
  
  para->geom->uniform = 0; //1: uniform; 0: non-uniform 
  para->mytime->dt = 0.1f; 
  para->mytime->t_steady = 100.0f;
  para->mytime->t_output =1000 ;
  para->prob->nu    = 1.79e-5f; //1.53e-5f;//1.53e-5f;
  para->prob->tur_model = LAM; //LAM, CHEN, 
  para->prob->coeff_h=0.004;//0.004
  para->prob->chen_a = 0.03874f; //0.03874f;  /* the coeffcient of Chen's model*/
  para->solv->solver = GS;
  para->solv->check_residual = 0; 
  para->outp->version = DEBUG; //DEMO, DEBUG;
  para->prob->Temp_opt=295;        //295.0f;
  para->prob->Prt=0.85f;
  para->bc->qs=0;
  para->prob->plume_mod=1;
  para->geom->plmax=3; //5
  para->geom->zv=-1.0f;

} // End of input_para( )