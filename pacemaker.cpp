#include "EDM_Struct.h"

/*******************************************************************************
                    STANDARD PARAMETERS FOR GENERAL MOTION
*******************************************************************************/
const double w = 2.0*M_PI;
const double q = 0.0;
const double F = 0.0;

const int dim = 2;
const int STEPS = 1000000;
const double dt = 0.0001;
const double t0 = 0.0;
const double x0 = 1.0;
const double v0 = 0.0;
/*******************************************************************************
                      AUXILIAR FUNCTIONS: DECLARATIONS
*******************************************************************************/
double steady_frec(EOM_Struct system, double f_top);

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
  pacemaker.Force.w = 0.15*w;
  pacemaker.Force.q = 1.5*q;
/*******************************************************************************
                         SOLVING EOM FOR GIVEN FORCE
*******************************************************************************/
  pacemaker.num_solve();
  pacemaker.dft_spectra(0.0);
  //pacemaker.DynSys.print_motion(0.0);

  double f_top = steady_frec(pacemaker, 5.0);
  printf("Force: %4.4f\tMain Frec: %4.4f\n",pacemaker.Force.F,f_top);
  /*
  std::cout << pacemaker.Force.w/w  << "\t"
            << pacemaker.Force.q  << "\t"
            << pacemaker.DynSys.dt  << "\t"
            << pacemaker.DynSys.NSTEP  << "\n";
  */
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
