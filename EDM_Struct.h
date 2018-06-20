#include "integration.h"
#include <cmath>
#include <fftw3.h>

/******************************************************************************/
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
                        STRUCTURES FOR EOM SOLVING
*******************************************************************************/

struct EOM_Force{
  double F;
  double w;
  double q;

  double force(int comp, double t, double *y);
};

struct EOM_Struct{
  EOM_Data DynSys;
  EOM_Force Force;

  void num_solve();
};

/*******************************************************************************
                      DEFINITION OF EXTERNAL FORCES
*******************************************************************************/

double EOM_Force::force(int comp, double t, double *y){
  if(comp == 1) return -w*w*y[0]-q*y[1]+F;
    else if (comp == 0) return y[1];
  else{
    std::cerr << "No more dependent variables" << '\n';
    return 0.0;
  }
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
