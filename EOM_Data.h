#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>
/*******************************************************************************
- AUTHOR: Diego Alejandro Herrera Rojas
- DATE: 07/08/18
- DESCRIPTION: Structure numerical integration, for project of
               Chaos Group (National University of Colombia). Remember to
               include this file with your routine-implementation files.
*******************************************************************************/
/*******************************************************************************
                       BASIC EOM DATA STORING STRUCTURE
*******************************************************************************/
struct EOM_Data{

  int dim;
  int NSTEP;
  double dt;
  double t0;
  double *init_data;
  double *motion;
  double *pdf;

  EOM_Data() {
    dim = 2;
    NSTEP = 10000;
    dt = 0.001;
    t0 = 0.0;
  };

  void initialize(int n, int N, double delta, double t_init);
  void kill();
  void WriteCoord(int i, int j, double x);
  double ReadCoord(int i, int j);
  void print_motion(double t_begin);
  void print_spectra(double f_top);

};
/*******************************************************************************
                       BASIC EOM DATA STORING STRUCTURE
*******************************************************************************/
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
  //free(pdf);
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
      printf("%4.7f\t %4.7e\t ",t,theta);
      for(int jj = 1; jj < dim; jj++)
        printf("%4.7e\t ",EOM_Data::ReadCoord(ii,jj));
      printf("\n");
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
