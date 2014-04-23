///////////////////////////////////////////////////////////////////////////////
///
/// \file   data_structure.h
///
/// \brief  Defining the the data structure and macro used in FFD
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
/// This file provides data structues used for storing problem properties, fluid
/// properties, simulation properties.
///
///////////////////////////////////////////////////////////////////////////////


#include <time.h>

#define IX(i,j,k) ((i)+(IMAX)*(j)+(IJMAX)*(k))+(IJK)
#define FIX(i,j,k) ((i)+(IMAX)*(j)+(IJMAX)*(k))
#define PIX(i,j,k) ((i)+ 2*(j)+4*(k))
#define FOR_EACH_CELL for(i=1; i<=imax; i++) { for(j=1; j<=jmax; j++) { for(k=1; k<=kmax; k++) {
#define FOR_ALL_CELL for(k=0; k<=kmax+1; k++) { for(j=0; j<=jmax+1; j++) { for(i=0; i<=imax+1; i++) {
#define FOR_U_CELL for(k=1; k<=kmax; k++) { for(j=1; j<=jmax; j++) { for(i=1; i<=imax-1; i++) {
#define FOR_V_CELL for(i=1; i<=imax; i++) { for(j=1; j<=jmax-1; j++) { for(k=1; k<=kmax; k++) {
#define FOR_W_CELL for(i=1; i<=imax; i++) { for(j=1; j<=jmax; j++) { for(k=1; k<=kmax-1; k++) {

#define FOR_KI for(i=1; i<=imax; i++) { for(k=1; k<=kmax; k++) {{
#define FOR_IJ for(i=1; i<=imax; i++) { for(j=1; j<=jmax; j++) {{
#define FOR_JK for(j=1; j<=jmax; j++) { for(k=1; k<=kmax; k++) {{
#define END_FOR }}}

#define REAL float
#ifndef max
  #define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#endif

#ifndef min
  #define min( a, b ) ( ((a) < (b)) ? (a) : (b) )
#endif


#define PI 3.1415926

#define X     0
#define Y     1
#define Z     2
#define VX    3
#define VY    4
#define VZ    5
#define VXM   6
#define VYM   7
#define VZM   8
#define VXS   9
#define VYS   10
#define VZS   11
#define DEN   12
#define DENS  13
#define IP    14
#define TMP1  15
#define TMP2  16
#define TMP3  17
#define TEMP  18
#define TEMPS 19
#define TEMPM 20
#define AP    21
#define AN    22
#define AS    23
#define AW    24
#define AE    25
#define AF    26
#define AB    27
#define  B    28
#define GX    29
#define GY    30
#define GZ    31
#define AP0   32
#define PP    33
#define FLAGP 34
#define FLAGU 35
#define FLAGV 36
#define FLAGW 37
#define LOCMIN 38
#define LOCMAX 39
#define VXBC 40
#define VYBC 41
#define VZBC 42
#define TEMPBC 43
#define FX 44
#define FY 45
#define FZ 46
#define ER 47
#define RR 48
#define VXC 49
#define VYC 50
#define VZC 51
#define FLUX 52
#define VSX 53
#define VSY 54
#define VSZ 55
#define CFLX 56
#define CFLY 57
#define CFLZ 58
#define DIST 59
#define VT   60
#define FLAGT   61


typedef enum{NOSLIP, SLIP, INFLOW, OUTFLOW, PERIODIC,SYMMETRY} BCTYPE;

typedef enum{TCONST, QCONST, ADIBATIC} BCTTYPE;

typedef enum{GS, TDMA} SOLVERTYPE;

typedef enum{SEMI, LAX, UPWIND, UPWIND_NEW} ADVECTION;

typedef enum{LAM, CHEN, CONST} TUR_MODEL;

typedef enum{BILINEAR,HYBRID, FSJ} INTERPOLATION;

typedef enum{DEMO, DEBUG} VERSION;

typedef enum{GRAVX,GRAVXN,GRAVY,GRAVYN,GRAVZ,GRAVZN} GRAVDIR;

typedef struct 
{
  REAL  Lx;       ///< domain size in x-direction                             
  REAL  Ly;       ///< domain size in y-direction                            
  REAL  Lz;       ///< domain size in z-direction                             
  int   imax;     ///< number of interior cells in x-direction               
  int   jmax;     ///< number of interior cells in y-direction               
  int   kmax;     ///< number of interior cells in z-direction               
  int   index[10];///< index of different boundary type
  int   zone_num; ///< total number of boundary cells
  int   zone_inlet;///< index number of inlet
  int   zone_outlet;///< index number of outlet
  int   zone_bl;///< index number of blockage
  int   zone_wall;///< index number of wall
  int   zone_us;///< index number of user defined sources
  int   zone_s; ///< index number of openings
  int   iplume; ///< I-index number of plume cell
  int   jplume;///< J-index number of plume cell
  int   kplume;///< K-index number of plume cell
  REAL  zv;    ///< virtual start location of plume model
  REAL  dx;    ///< length delta_x of one cell in x-direction            
  REAL  dy;    ///< length delta_y of one cell in y-direction              
  REAL  dz;    ///<length delta_z of one cell in z-direction              
} GEOM_DATA;

typedef struct{
  int   cal_mean;    ///< 1: calculate mean value                              
  REAL  v_ref;       ///< reference velocity                                   
  REAL  Temp_ref;    ///< reference temperature                                
  int   plot_grid;   ///< number of plotting grids                               
  REAL  v_length;    ///< the ratio factor of the velocity length              
  int   i_N;         ///< the number of grids plotted in x direction           
  int   j_N;         ///< the number of grids plotted in y direction           
  int   winx;        ///< the resolution of screen at x direction              
  int   winy;       ///< the resolution of screen at y direction             
  VERSION  version;  ///< DEMO, DEBUG                                              
} OUTP_DATA;

typedef struct{
  REAL  RE;       ///< Reynolds number Re                                     
  REAL  mu;       ///< physical viscosity                                     
  REAL  nu;       ///<  kinematic viscosity                                    
  REAL  rho;      ///<  density                                                
  REAL  diff;     ///<  diffusivity for particle density                       
  REAL  alpha;    ///<  thermal diffusity                                      
  REAL  alpha_co; ///<  proption coefficient of thermal diffusity at the wall  
  REAL  coeff_h;  ///<  convective coefficient
  REAL  k;        ///<  thermal conductivity                                   
  REAL  gravx;    ///< gravity coefficient in x direction
  REAL  gravy;     ///< gravity coefficient in y direction
  REAL  gravz;     ///< gravity coefficient in z direction   
  GRAVDIR gravdir;  ///< gravity direction
  REAL  beta;     ///<  coefficient of thermal expansion                       
  REAL  cond;     ///<  conductivity 
  REAL  trefmax;  ///< reference temperature 
  REAL  spec;     ///< specific heat 
  int   readfile; ///< Read old data file as initial value(1:yes, 0:no)       */
  int   moive;    ///< output data for make animation file(1:yes, 0:no)       */
  int   output;   ///< 0: have not been written; 1: done                      */ 
  int   plume_mod; ///< 0: not use; 1: use                                    */ 
  TUR_MODEL tur_model; ///< LAM, CHEN, 100NU                                   */ 
  REAL  chen_a;   ///< coefficeint of Chen's zero euqation turbulence model   */
  REAL  Prt;      ///< turbulent Prandl number */
  REAL  Temp_opt; ///< operational temperature
  REAL  resu;     ///< residual of u 
  REAL  resv;      ///< residual of v 
  REAL  resw;      ///< residual of w 
  REAL  resp;      ///< residual of p 
  int  iteru;      ///< iteration number of u 
  int  iterv;      ///< iteration number of v 
  int  iterw;      ///< iteration number of w 
  int  iterp;      ///< iteration number of p 
}PROB_DATA;

typedef struct{
  int   NBOUT; ///< number of outflow boundary
  REAL  qs;    ///< heat generation rate of heat source cell
  REAL  u_bc[200];///< u boundary value
  REAL  v_bc[200];///< v boundary value
  REAL  w_bc[200];///< w boundary value
  REAL  t_bc[200];///< Temperature boundary value
  REAL  d_bc[200];///< concentration boundary value
  REAL  fltmp[200];///< Temperature boudnary condition type 
  REAL  um_bc[200];///< momentum source value at the bounary 
  REAL  zv_bc[200];///< virtual location value of heat source cell
}BC_DATA;

typedef struct 
{
  REAL   dt;         ///<  time step size                                      
  REAL   t;          ///<  current time                                        
  REAL   t_steady;   ///< necessary time for steady flow                      
  int     t_output;  ///<  the interval of iteration step to output data       
  int     t_step;     ///<  current iteration step                              
  clock_t t_start;    ///<  starting CPU time                                   
  clock_t t_end;      ///<  ending CPU time                                     
}TIME_DATA;

typedef struct 
{
  int   caseID;       ///<  1: Pure Conduction with uniform grid                
  REAL   f_dif;      ///<  f_dif=0: explict scheme; f_dif=1: implict scheme    
  SOLVERTYPE solver;  ///<  GS, TDMA                                            
  int check_residual; ///<  1: check, 0: donot check                            
  ADVECTION advection_solver; ///<  SEMI, LAX, UPWIND, UPWIND_NEW */  
  INTERPOLATION interpolation; ///<  BILINEAR, FSJ */
  int  read_file;     ///<  1: Read previous file                              
}SOLV_DATA;

typedef struct 
{
  GEOM_DATA  *geom; ///< Pointer to geometry parameter
  OUTP_DATA  *outp; ///< Pointer to output parameter
  PROB_DATA  *prob; ///< Pointer to problem parameter
  TIME_DATA  *mytime;///< Pointer to time parameter
  BC_DATA    *bc;    ///< Pointer to boundary parameter
  SOLV_DATA  *solv;  ///< Pointer to solver parameter
}PARA_DATA;
