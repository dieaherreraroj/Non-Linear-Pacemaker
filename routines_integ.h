#include "integration.h"
#include <cmath>
#include <fftw3.h>
//#include "omp.h"

/*******************************************************************************
                        AUXILIAR FUNCTIONS: DECLARATIONS
*******************************************************************************/

void copy_vect(int n, double *x, double *y);
void update_kvect(int n, int tag, fptr f, double h, double t, double *y, double *ki, double *kf);
void aux_interm(int n, int tag, double h, double *x, double *y, double *z);

/*******************************************************************************
                     INTEGRATION ALGORITHM: RUNGE-KUTTA(4)
*******************************************************************************/

void EOM_Data::rk4_integration(fptr f){
  // Allocate memory
  double *y = (double*) calloc(dim,sizeof(double));
  double *k1 = (double*) calloc(dim,sizeof(double));
  double *k2 = (double*) calloc(dim,sizeof(double));
  double *k3 = (double*) calloc(dim,sizeof(double));
  double *k4 = (double*) calloc(dim,sizeof(double));
  double t = t0;
  // Set initial conditions
  copy_vect(dim,init_data,y);
  // Core of routine: Integration steps
  for(int ii = 0; ii < NSTEP; ii++){
    t = t0 + ii*dt;
    // Update numerical solution data in EOM_Data Structure
    WriteCoord(ii,0,atan2(sin(y[0]),cos(y[0])));
    for(int jj = 1; jj < dim; jj++)
      WriteCoord(ii,jj,y[jj]);
    // Carry out next integration step
    update_kvect(dim,1,f,dt,t,y,y,k1);
    update_kvect(dim,2,f,dt,t,y,k1,k2);
    update_kvect(dim,3,f,dt,t,y,k2,k3);
    update_kvect(dim,4,f,dt,t,y,k3,k4);
    for(int jj = 0; jj<dim; jj++)
      y[jj] += (1.0/6.0)*dt*(k1[jj] + k4[jj] + 2.0*(k2[jj] + k3[jj]));
  }
  // Release memory
  free(y);
  free(k1);
  free(k2);
  free(k3);
  free(k4);
}

/*******************************************************************************
                       AUXILIAR FUNCTIONS: IMPLEMENTATIONS
*******************************************************************************/

void copy_vect(int n, double *x, double *y){
  for(int ii = 0; ii<n; ii++)
    y[ii] = x[ii];
}

/******************************************************************************/

void aux_interm(int n, int tag, double h, double *x, double *y, double *z){
  if(tag == 1)
    copy_vect(n,y,z);
  if(tag == 2 || tag == 3)
    for(int ii = 0; ii<n; ii++)
      z[ii] = x[ii] + 0.5*h*y[ii];
  if(tag == 4)
    for(int ii = 0; ii<n; ii++)
      z[ii] = x[ii] + h*y[ii];
}

/******************************************************************************/

void update_kvect(int n, int tag, fptr f, double h, double t, double *y, double *ki, double *kf){
  double *aux = (double*) calloc(n,sizeof(double));
  double t_aux = 0.0;
  if(tag == 1) t_aux = t;
    else if(tag == 2 || tag == 3) t_aux = t + 0.5*h;
      else if(tag == 4) t_aux = t + h;
  else t_aux = 0.0;
  aux_interm(n,tag,h,y,ki,aux);
  for(int ii = 0; ii<n; ii++)
    kf[ii] = f(ii,t_aux,aux);
  free(aux);
}
