#include "integration.h"

/*******************************************************************************
Please do not modify any of the routines below this barrier, otherwise, you may
have bad results because of bad implementations, only access the routine above,
or the file integration.h.

P. S. If it is necessary, check with Diego Alejandro Herrera Rojas (UNAL) via
      e-mail: dieaherreraroj@unal.edu.co
*******************************************************************************/
/******************************************************************************
                     BASIC EOM DATA STRUCTURE FUNCTIONS
******************************************************************************/

void EOM_Data::initialize(int n, int N, double delta, double t_init){
  if(0 < n && 0 < N && 0 < delta){
    NSTEP = N;
    dim = n;
    dt = delta;
    t0 = t_init;
    init_data = (double*) calloc(dim,sizeof(double));
    motion = (double*) calloc(dim*NSTEP,sizeof(double));
  }
  else
    std::cerr << "Inappropriate initialization" << '\n';
}

void EOM_Data::kill(){
  free(init_data);
  free(motion);
  free(pdf);
}

void EOM_Data::WriteCoord(int i, int j, double x){
  if(i < NSTEP && j < dim)
    motion[dim*i+j] = x;
  else
    std::cerr << "Writing Outside Array" << '\n';
}

double EOM_Data::ReadCoord(int i, int j){
  if(i < NSTEP && j < dim)
    return motion[dim*i+j];
  else{
    std::cerr << "Reading Outside Array" << '\n';
    return 0.0;
  }
}

void EOM_Data::print_motion(double t_begin){
  double t = 0.0;
  for(int ii = 0; ii < NSTEP; ii++){
    t = t0 + ii*dt;
    if(t > t_begin){
      double ang = EOM_Data::ReadCoord(ii,0);
      double theta = atan2(sin(ang),cos(ang));
      printf("%4.7f\t %4.7e\t %4.7e\n",t,theta,EOM_Data::ReadCoord(ii,1));
    }
  }
}

void EOM_Data::print_spectra(double f_top){
  double f = 0.0;
  // Strange feature. Requires further study
  double f_samp = 1.0/dt;
  for(int ii = 0; f < f_top && ii < NSTEP; ii++){
    f = (f_samp*ii)/NSTEP;
    printf("%4.7f\t %4.7e\n",f,pdf[ii]);
  }
}

/*******************************************************************************
               AUXILIAR FUNCTIONS (DECLARATION & IMPLEMENTATION)
*******************************************************************************/
void copy_vect(int n, double *x, double *y);
void aux_interm(int n, int tag, double h, double *x, double *y, double *z);
/******************************************************************************/

void copy_vect(int n, double *x, double *y){
  for(int ii = 0; ii<n; ii++)
    y[ii] = x[ii];
}

void aux_interm(int n, int tag, double h, double *x, double *y, double *z){
  if(tag == 2 || tag == 3)
    for(int ii = 0; ii<n; ii++)
      z[ii] = x[ii] + 0.5*h*y[ii];
  if(tag == 4)
    for(int ii = 0; ii<n; ii++)
      z[ii] = x[ii] + h*y[ii];
}

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
/*
void EOM_Struct::dft_spectra(double t_min){
  fftw_complex *in, *out;
  fftw_plan plan;
  // Allocate memory
  in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*DynSys.NSTEP);
  out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*DynSys.NSTEP);
  plan = fftw_plan_dft_1d(DynSys.NSTEP,in,out,FFTW_FORWARD,FFTW_ESTIMATE);
  // Fill with stationary state
  double t = DynSys.t0;
  for(int ii = 0; ii < DynSys.NSTEP; ii++){
    t = DynSys.t0 + ii*DynSys.dt;
    if(t_min <= t)
      in[ii][0] = DynSys.motion[ii];
    else
      in[ii][0] = 0.0;
    in[ii][1] = 0.0;
  }
  // Execute fft algorithm
  fftw_execute(plan);
  // Fill with estimated pdf
  DynSys.pdf = (double*) calloc(DynSys.NSTEP,sizeof(double));
  double p = 0.0;
  for(int ii = 0; ii < DynSys.NSTEP; ii++){
    p = out[ii][0]*out[ii][0] + out[ii][1]*out[ii][1];
    DynSys.pdf[ii] = DynSys.dt*DynSys.dt*p;
  }
  // Deallocate Memory
  fftw_destroy_plan(plan);
  fftw_free(in);
  fftw_free(out);
}
*/
