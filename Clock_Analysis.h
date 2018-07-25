#include "Clock_Struct.h"
//#include <fftw3.h>

/*******************************************************************************
CREATOR: Diego Alejandro Herrera Rojas.
FOR: Grupo Caos y Complejidad - National University of Colombia.
********************************************************************************
PROJECT: Clock Escapements as Non-linear Pacemakers.
DESCRIPTION: Implementation of functions used in numerical solution of equations
             of motion.
NOTES: Compile in UNIX system, with a compiler that supports C++11 (at least).
       Library fftw3 (last version), must be available.

       Compile using g++ -std=c++11 pacemaker.cpp -lfftw3
*******************************************************************************/
/*******************************************************************************
Please do not modify any of the routines below this barrier, otherwise, you may
have bad results because of bad implementations, only access the routine above,
or the file integration.h.

P. S. If it is necessary, check with Diego Alejandro Herrera Rojas (UNAL) via
      e-mail: dieaherreraroj@unal.edu.co
*******************************************************************************/
/*******************************************************************************
                        NUMERICAL INTEGRATION WITH rk4
*******************************************************************************/
void Clock::motion_integration(){
  // Allocate memory
  double *y = (double*) calloc(integ_data.dim,sizeof(double));
  double *k1 = (double*) calloc(integ_data.dim,sizeof(double));
  double *k2 = (double*) calloc(integ_data.dim,sizeof(double));
  double *k3 = (double*) calloc(integ_data.dim,sizeof(double));
  double *k4 = (double*) calloc(integ_data.dim,sizeof(double));
  double *aux = (double*) calloc(integ_data.dim,sizeof(double));
  double t = integ_data.t0;
  double phi = 2.0*M_PI/num_teeth;
  int new_tag = 0;
  int change_phase = 0;
  // Set initial conditions
  copy_vect(integ_data.dim,integ_data.init_data,y);
  double phi_n = floor(y[2]/phi)*phi;
  // Core of routine: Integration steps
  for(int ii = 0; ii < integ_data.NSTEP; ii++){
    t = integ_data.t0 + ii*integ_data.dt;

    // Record motion of clock : include constraints
    if(change_phase == 1 && new_tag%2 == 1){
      phi_n = floor(y[2]/phi)*phi;
    }
    if(new_tag%2 == 0){                       // uncoupled trnasition
      integ_data.EOM_Data::WriteCoord(ii,0,atan2(sin(y[0]),cos(y[0])));
      for(int jj = 1; jj < integ_data.dim; jj++)
        integ_data.EOM_Data::WriteCoord(ii,jj,y[jj]);
    }
    if(new_tag == 1){                         // coupled phase : constraints
      integ_data.EOM_Data::WriteCoord(ii,0,atan2(sin(y[0]),cos(y[0])));
      integ_data.EOM_Data::WriteCoord(ii,1,y[1]);
      integ_data.EOM_Data::WriteCoord(ii,2,phi_n + rho*(y[0]/tan(alpha)));
      integ_data.EOM_Data::WriteCoord(ii,1,y[1]);
    }
    if(new_tag == 3){                         // coupled phase : constraints
      integ_data.EOM_Data::WriteCoord(ii,0,atan2(sin(y[0]),cos(y[0])));
      integ_data.EOM_Data::WriteCoord(ii,1,y[1]);
      integ_data.EOM_Data::WriteCoord(ii,2,phi_n - rho*(y[0]/tan(alpha)));
      integ_data.EOM_Data::WriteCoord(ii,1,y[1]);
    }

    // Carry out integration step
    // Update k1
    for(int hh = 0; hh < integ_data.dim; hh++)
      k1[hh] = Clock::force(hh,t,y);
    // Update k2
    aux_interm(integ_data.dim,2,integ_data.dt,y,k1,aux);
    for(int hh = 0; hh < integ_data.dim; hh++)
      k2[hh] = Clock::force(hh,t+0.5*integ_data.dt,aux);
    // Update k3
    aux_interm(integ_data.dim,3,integ_data.dt,y,k2,aux);
    for(int hh = 0; hh < integ_data.dim; hh++)
      k3[hh] = Clock::force(hh,t+0.5*integ_data.dt,aux);
    // Update k4
    aux_interm(integ_data.dim,4,integ_data.dt,y,k3,aux);
    for(int hh = 0; hh < integ_data.dim; hh++)
      k4[hh] = Clock::force(hh,t+integ_data.dt,aux);
    // Fill array with data generated
    for(int jj = 0; jj < integ_data.dim; jj++)
      y[jj] += (1.0/6.0)*integ_data.dt*(k1[jj] + k4[jj] + 2.0*(k2[jj] + k3[jj]));

    // Check if integration step leads to change of phase
    double d = y[2] - floor(y[2]/phi);
    if(fabs(y[0]) <= (phi_max/tan(alpha))*rho){
      if(0 <= d && d <= phi_max){
        new_tag = 1;                            // First coupled phase (A1)
      }
      else if(phi_max < d && d <= phi/2.0){
        new_tag = 2;                            // First uncoupled phase (B1)
      }
      else if(phi/2.0 < d && d <= (phi_max + phi/2.0)){
        new_tag = 3;                            // Second coupled phase (A2)
      }
      else if((phi_max + phi/2.0) < d && d <= phi){
        new_tag = 4;                            // Second uncoupled phase (B2)
      }
      else{
        std::cerr << "Outside allowed range of Cogwheel." << '\n';
        break;
      }
    }
    else{
      std::cerr << "Outside allowed range of Pendullum." << '\n';
      break;
    }

    // Emulate collision cogwheel - pendullum : impose boundary conditions
    if(tag_phase != new_tag){
      change_phase = 1;
      if(new_tag == 1 || new_tag == 3){
        y[3] = y[1]*tan(alpha);
      }
      else{
        continue;
      }
    }
    else{
      change_phase = 0;
      continue;
    }
  }
  // Release memory
  free(y);
  free(k1);
  free(k2);
  free(k3);
  free(k4);
  free(aux);
}
