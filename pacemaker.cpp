#include "EDM_Struct.h"

const double w = 2.0*M_PI;
const double q = 0.0;
const double F = 0.0;

const int dim = 2;
const int STEPS = 10000;
const double dt = 0.001;
const double t0 = 0.0;
const double x0 = 1.0;
const double v0 = 0.0;

int main(void){

  EOM_Struct pacemaker;

/*******************************************************************************
                         SETTING UP FORCE PARAMETERS
*******************************************************************************/
  pacemaker.Force.F = 2.5*F;
  pacemaker.Force.w = 1.0*w;
  pacemaker.Force.q = 1.5*q;
/*******************************************************************************
                        SETTING UP INITIAL CONDITIONS
*******************************************************************************/
  pacemaker.DynSys.initialize(dim,STEPS,dt,t0);
  pacemaker.DynSys.init_data[0] = x0;
  pacemaker.DynSys.init_data[1] = v0;
/*******************************************************************************
                         SOLVING EOM FOR GIVEN FORCE
*******************************************************************************/
  pacemaker.num_solve();
  pacemaker.DynSys.print_motion(t0);
/******************************************************************************/
  pacemaker.DynSys.kill();
}
