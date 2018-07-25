#include "EOM_Struct.h"

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
             ANALISYS ROUTINES: EOM SOLVING AND FOURIER ANALISYS
********************************************************************************
                            EDM NUMERICAL SOLVING
*******************************************************************************/

void EOM_Struct::num_solve(){
  // Allocate memory
  double *y = (double*) calloc(DynSys.dim,sizeof(double));
  double *k1 = (double*) calloc(DynSys.dim,sizeof(double));
  double *k2 = (double*) calloc(DynSys.dim,sizeof(double));
  double *k3 = (double*) calloc(DynSys.dim,sizeof(double));
  double *k4 = (double*) calloc(DynSys.dim,sizeof(double));
  double *aux = (double*) calloc(DynSys.dim,sizeof(double));
  double t = DynSys.t0;
  // Set initial conditions
  copy_vect(DynSys.dim,DynSys.init_data,y);
  // Core of routine: Integration steps
  for(int ii = 0; ii < DynSys.NSTEP; ii++){
    t = DynSys.t0 + ii*DynSys.dt;
    // Update numerical solution data in EOM_Data Structure
    DynSys.EOM_Data::WriteCoord(ii,0,atan2(sin(y[0]),cos(y[0])));
    for(int jj = 1; jj < DynSys.dim; jj++)
      DynSys.EOM_Data::WriteCoord(ii,jj,y[jj]);
    // Carry out next integration step
    // Update k1
    for(int hh = 0; hh < DynSys.dim; hh++)
      k1[hh] = Force.EOM_Force::force(hh,t,y);
    // Update k2
    aux_interm(DynSys.dim,2,DynSys.dt,y,k1,aux);
    for(int hh = 0; hh < DynSys.dim; hh++)
      k2[hh] = Force.EOM_Force::force(hh,t+0.5*DynSys.dt,aux);
    // Update k3
    aux_interm(DynSys.dim,3,DynSys.dt,y,k2,aux);
    for(int hh = 0; hh < DynSys.dim; hh++)
      k3[hh] = Force.EOM_Force::force(hh,t+0.5*DynSys.dt,aux);
    // Update k4
    aux_interm(DynSys.dim,4,DynSys.dt,y,k3,aux);
    for(int hh = 0; hh < DynSys.dim; hh++)
      k4[hh] = Force.EOM_Force::force(hh,t+DynSys.dt,aux);
    // Fill array with data generated
    for(int jj = 0; jj < DynSys.dim; jj++)
      y[jj] += (1.0/6.0)*DynSys.dt*(k1[jj] + k4[jj] + 2.0*(k2[jj] + k3[jj]));
  }
  // Release memory
  free(y);
  free(k1);
  free(k2);
  free(k3);
  free(k4);
  free(aux);
}

/*******************************************************************************
                      STEADY STATE MOTION FOURIER ANALISYS
*******************************************************************************/

void EOM_Struct::dft_spectra(double t_min){
  // Variables for Fourier Analysis
  fftw_complex *in, *out;
  fftw_plan plan;
  // Allocate memory
  DynSys.pdf = (double*) calloc(DynSys.NSTEP, sizeof(double));
  in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*DynSys.NSTEP);
  out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*DynSys.NSTEP);
  // Create plan
  plan = fftw_plan_dft_1d(DynSys.NSTEP,in,out,FFTW_FORWARD,FFTW_ESTIMATE);
  // Test Filling: Works ok. Problem is in array filling.
  double t = 0.0;
  for(int ii = 0; ii < DynSys.NSTEP; ii++){
    t = DynSys.t0 + ii*DynSys.dt;
    /*
    in[ii][0] = cos(Force.w*t);
    in[ii][1] = 0.0;*/

    if(t >= t_min)
      in[ii][0] = DynSys.EOM_Data::ReadCoord(ii,0);
    else
      in[ii][0] = 0.0;
    in[ii][1] = 0.0;
  }
  // Compute DFT
  fftw_execute(plan);
  double p = 0.0;
  // Fill dft array
  for(int ii = 0; ii < DynSys.NSTEP; ii++){
    p = out[ii][0]*out[ii][0] + out[ii][1]*out[ii][1];
    DynSys.pdf[ii] = DynSys.dt*DynSys.dt*p;
  }
}
