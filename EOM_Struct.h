#include "EOM_Data.h"
#include <fftw3.h>

/*******************************************************************************
Please, just modify parameters of structure EOM_Force if you need to solve a
different equation of motion, with different parameters. Go to end of file if
you want to change the Force.

P. S. Mantaining the EDM_Struct and EOM_Data struct is essencial for the good
      behavior of the programm. Under no circumstances should you modify these
      structures.

In this file you have acces to all declarations of functions of importance. They
tell you wich are the parameters of each one.
*******************************************************************************/
/*******************************************************************************
                             AUXILIAR DEFINITIONS
*******************************************************************************/

using fptr = double(int, double, double *);

/*******************************************************************************
                        STRUCTURES FOR EOM SOLVING
********************************************************************************
This is the only part you should modify, in here you can define all parameters
for almos every external force you can imagine. Please do not try to change any
other part of the file. Define the parameters inside EOM_Force as desired and
remember to modify memeber function force() accordingly.
*******************************************************************************/

struct EOM_Force{
  double a;
  double It;
  double Ip;
  double gamma;
  double alfa;
  double Itot;
  double w;

  EOM_Force() {
    a = 1.0;
    It = 1.0;
    Ip = 0.1;
    gamma = 0.2;
    alfa = tan(14.0*M_PI/180.0);
    Itot = It + Ip/(alfa*alfa);
    w = 2.0*M_PI;
  };

  double force(int comp, double t, double *y);
};

/******************************************************************************/

struct EOM_Struct{
  EOM_Data DynSys;
  EOM_Force Force;

  void num_solve();
  void dft_spectra(double t_min);
};

/*******************************************************************************
                      DEFINITION OF EXTERNAL FORCES
********************************************************************************
This you can modify as much as you desire, this function contains information of
external forces on the system. For sake of consistency, you can modify force pa-
rameters by editing EOM_Force structure.
*******************************************************************************/

double EOM_Force::force(int comp, double t, double *y){
  if(comp == 0) return y[1];
    else if (comp == 1) return -(gamma/Itot)*y[1]-(It/Itot)*w*w*y[0]+a/alfa;
  else{
    std::cerr << "No more dependent variables" << '\n';
    return 0.0;
  }
}
