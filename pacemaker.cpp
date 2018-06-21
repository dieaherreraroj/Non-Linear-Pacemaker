#include "EDM_Struct.h"

/*******************************************************************************
                    STANDARD PARAMETERS FOR GENERAL MOTION
*******************************************************************************/
const double w = 2.0*M_PI;
const double q = 0.0;
const double F = 0.0;

const int dim = 2;
const int STEPS = 100000;
const double dt = 0.001;
const double t0 = 0.0;
const double x0 = 1.0;
const double v0 = 0.0;
/*******************************************************************************
                      AUXILIAR FUNCTIONS: DECLARATIONS
*******************************************************************************/
double steady_frec(EOM_Struct system);

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
  pacemaker.Force.F = 2.5*F;
  pacemaker.Force.w = 1.0*w;
  pacemaker.Force.q = 1.5*q;
/*******************************************************************************
                         SOLVING EOM FOR GIVEN FORCE
*******************************************************************************/
  pacemaker.num_solve();
  pacemaker.dft_spectra(0.0);
  pacemaker.DynSys.print_spectra(5.0);
  double f_top = steady_frec(pacemaker);
  //printf("Force: %4.4f\tMain Frec: %4.4f\n",pacemaker.Force.F,f_top);
/*******************************************************************************
                      RETURN MEMORY (DO NOT DELETE !!!)
*******************************************************************************/
  pacemaker.DynSys.kill();
}

/*******************************************************************************
                     AUXILIAR FUNCTIONS: IMPLEMENTATIONS
*******************************************************************************/
double steady_frec(EOM_Struct system){
  int peak = 0;
  for(int ii = 0; ii < system.DynSys.NSTEP; ii++){
    if(system.DynSys.pdf[peak] < system.DynSys.pdf[ii])
      peak = ii;
  }
  double f_samp = 2.0/system.DynSys.dt;
  return (f_samp*peak)/system.DynSys.NSTEP;
}
