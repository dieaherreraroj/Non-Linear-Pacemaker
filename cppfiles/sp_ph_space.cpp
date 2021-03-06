#include "EOM_Analysis.h"

/*******************************************************************************
CREATOR: Diego Alejandro Herrera Rojas.
FOR: Grupo Caos y Complejidad - National University of Colombia.
********************************************************************************
PROJECT: Clock Escapements as Non-linear Pacemakers.
DESCRIPTION: Numerical solution of equations of motion for different dynamical
             systems, including a Fourier Analysis. Min function, used for par-
             ticular simlations and data generation.
NOTES: Compile in UNIX system, with a compiler that supports C++11 (at least).
       Library fftw3 (last version), must be available.

       Compile using g++ -std=c++11 pacemaker.cpp -lfftw3
*******************************************************************************/
/*******************************************************************************
                    STANDARD PARAMETERS FOR GENERAL MOTION
*******************************************************************************/
const double w = 1.0;
const double q = 0.25;
const double F = 1.0;

const int dim = 2;
const int STEPS = 1000 * ((int) (100.0/q)) + 1000000;
const double dt = 0.001;
const double t0 = 0.0;
const double x0 = 1.0;
const double v0 = 0.0;
/*******************************************************************************
                      AUXILIAR FUNCTIONS: DECLARATIONS
*******************************************************************************/
double steady_frec(EOM_Struct system, double f_top);
/*******************************************************************************
                     MAIN FUNCTION: REGION FOR WORKING
*******************************************************************************/
int main(void){

  EOM_Struct pacemaker;
/*******************************************************************************
                        SETTING UP INITIAL CONDITIONS
*******************************************************************************/
  pacemaker.DynSys.initialize(dim,STEPS,dt,t0);
  pacemaker.DynSys.init_data[0] = x0;
  pacemaker.DynSys.init_data[1] = v0;
/*******************************************************************************
                         SETTING UP FORCE PARAMETERS
*******************************************************************************/
  pacemaker.Force.w = w;
  pacemaker.Force.q = q;
  pacemaker.Force.F = 1.5*F;
/******************************************************************************/
/*******************************************************************************
                         SOLVING EOM FOR GIVEN FORCE
*******************************************************************************/
  pacemaker.num_solve();
  pacemaker.DynSys.print_motion(100.0/q);
/*******************************************************************************
                      RETURN MEMORY (DO NOT DELETE !!!)
*******************************************************************************/
  pacemaker.DynSys.kill();
}

/*******************************************************************************
                     AUXILIAR FUNCTIONS: IMPLEMENTATIONS
*******************************************************************************/
double steady_frec(EOM_Struct system, double f_top){
  int peak = 0;
  double f_samp = 1.0/system.DynSys.dt;
  double f = 0.0;
  for(int ii = 0; ii < system.DynSys.NSTEP && f <= f_top; ii++){
    f = (f_samp*ii)/system.DynSys.NSTEP;
    if(system.DynSys.pdf[peak] < system.DynSys.pdf[ii])
      peak = ii;
  }
  if(0 <= peak)
    return (f_samp*peak)/system.DynSys.NSTEP;
  else
    return -1.0;
}
