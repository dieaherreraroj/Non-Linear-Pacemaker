#include "routines_integ.h"

const double w = 2.0*M_PI;
const double q = 0.0;
const double F = 0.0;

double force(int comp, double t, double *y);

int main(void){

  const int dim = 2;
  const int STEPS = 10000;
  const double dt = 0.001;
  const double t0 = 0.0;
  const double x0 = 1.0;
  const double v0 = 0.0;

  EOM_Data pacemaker;
  pacemaker.initialize(dim,STEPS,dt,t0);

  pacemaker.init_data[0] = x0;
  pacemaker.init_data[1] = v0;

  pacemaker.rk4_integration(force);
  double t = 0.0;

  for(int ii = 0; ii < pacemaker.NSTEP; ii++){
    t = pacemaker.t0 + ii*pacemaker.dt;
    //if(t > 70.0/q){
      double ang = pacemaker.ReadCoord(ii,0);
      double theta = atan2(sin(ang),cos(ang));
      printf("%4.7f\t %4.7e\t %4.7e\n",t,theta,pacemaker.ReadCoord(ii,1));
    //}
  }

}

double force(int comp, double t, double *y){
  if(comp == 1) return -w*w*y[0]-q*y[1]+F;
    else if (comp == 0) return y[1];
  else{
    std::cerr << "No more dependent variables" << '\n';
    return 0.0;
  }
}
