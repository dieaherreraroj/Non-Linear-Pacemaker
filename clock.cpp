#include "Clock_Analysis.h"

/*******************************************************************************
                       BASIC DYNAMICAL AND DESIGN CONSTANTS
*******************************************************************************/
const double q = 0.25;
const double w = 2.0*M_PI;
const double alpha = 1.39;
const double phi = 0.25;
const double I = 1.0;
const double tau = 20.0;
const double r = 1.0;
/*******************************************************************************
                         NUMERICAL INTEGRATION CONSTANTS
*******************************************************************************/
const int dim = 4;
const double t0 = 0.0;
const double tf = 100.0/q;
const double dt = 0.001;

int main(void){
  Clock grandpa;
  int NSTEP = (int) ((tf-t0)/dt);
/*******************************************************************************
                      SET CLOCK DESIGN & MECHANICAL PARAMS.
*******************************************************************************/
  grandpa.init_design_param(phi,r,alpha,15,q);
  grandpa.init_pend_mec_param(w,I,0.1*q);
  grandpa.init_wheel_mec_param(tau,0.25*I,0.1*q);
/*******************************************************************************
                           SET INTEGRATION PARAMETERS
*******************************************************************************/
  grandpa.integ_data.initialize(dim,NSTEP,dt,t0);
/*******************************************************************************
                           SET INITIAL CONDITIONS
*******************************************************************************/
  for(int ii = 0; ii < grandpa.integ_data.dim; ii++)
    grandpa.integ_data.init_data[ii] = 0.0;
/*******************************************************************************
                           SOLVE EQUATIONS OF MOTION
*******************************************************************************/
  grandpa.motion_integration();
  grandpa.integ_data.print_motion(t0);
/*******************************************************************************
                                 RETURN MEMORY
*******************************************************************************/
  free(grandpa.integ_data.init_data);
  free(grandpa.integ_data.motion);
}
