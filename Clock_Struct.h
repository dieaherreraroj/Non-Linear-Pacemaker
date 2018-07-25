#include "EOM_Data.h"
#include <fftw3.h>

/*******************************************************************************
CREATOR: Diego Alejandro Herrera Rojas.
FOR: Grupo Caos y Complejidad - National University of Colombia.
********************************************************************************
PROJECT: Clock Escapements as Non-linear Pacemakers.
DESCRIPTION: Model of clocks as two-dof systems, with a complicated external
             force, involving coupled-uncoupled regimes. Structure and memeber
             function declarations.
NOTES: Compile in UNIX system, with a compiler that supports C++11 (at least).
       Library fftw3 (last version), must be available.
       Remember to include this header in your program.
       Compile using g++ -std=c++11 pacemaker.cpp -lfftw3
*******************************************************************************/
/*******************************************************************************
                                CLOCK STRUCTURE
*******************************************************************************/
struct Clock{
  int tag_phase;                      // Check if coupled/uncoupled motion
/*******************************************************************************
                               DESIGN PARAMETERS
*******************************************************************************/
  double phi_max;                     // Maximum angle of coupled motion
  double rho;                         // Ratio pendullum - cogwheel radii
  double alpha;                       // Slope parameter
  int num_teeth;                      // Number of cogwheel teeth
  double int_fric;                    // Internal friction while coupled motion
/*******************************************************************************
                             MECHANICAL PARAMETERS
*******************************************************************************/
  // Pendullum parameters
  double pend_fric;                   // Friction of oscillatory motion
  double pend_inertia;                // Moment of inertia of pendullum
  double omega;                       // Natural frequency of oscillation
  // Cogwheel parameters
  double wheel_fric;                  // Friction of cogwheel rotation
  double wheel_inertia;               // Moment of inertia of cogwheel
  double torque;                      // Constant external torque
/*******************************************************************************
                             DEFAULT CONSTRUCTOR
*******************************************************************************/
  Clock(){
    rho = 1.0;
    alpha = 0.1;
    num_teeth = 15;
    int_fric = 0.25; pend_fric = 0.25; wheel_fric = 0.25;
    pend_inertia = 4.17e-5; wheel_inertia = 8.34e-5;
    omega = 24.248;
    torque = 10.0;
  }
/*******************************************************************************
                          INTEGRATION DATA AND PARAM.
*******************************************************************************/
  EOM_Data integ_data;
/*******************************************************************************
                         INITIALIZING MEMBER FUNCTIONS
*******************************************************************************/
  void init_design_param(double phi, double r, double a, double teeth, double fric);
  void init_pend_mec_param(double w, double I, double fric);
  void init_wheel_mec_param(double tau, double I, double fric);
/*******************************************************************************
                          EXTERNAL FORCES OF SYSTEM
*******************************************************************************/
  double force (int comp, double t, double *y);
/*******************************************************************************
                     MOTION INTEGRATION MEMBER FUNCTION
*******************************************************************************/
  void motion_integration();
};
/*******************************************************************************
                          EXTERNAL FORCES OF SYSTEM
*******************************************************************************/
double Clock::force(int comp, double t, double *y){
  double phi = 2.0*M_PI/num_teeth;
  double d = y[2] - floor(y[2]/phi);
  if(fabs(y[0]) <= (phi_max/tan(alpha))*rho){
    if(0 <= d && d <= phi_max){
      tag_phase = 1;                            // First coupled phase (A1)
      double tang = tan(alpha);
      double tang2 = tang*tang;
      double I = (pend_inertia*tang2 + wheel_inertia)/tang2;
      double Tr = torque*(wheel_inertia/(I*tang));
      double Io = pend_inertia/I;
      double Ir = wheel_inertia/I;
      if (comp == 0) return y[1];
      else if (comp == 1) return -Io*omega*omega*y[0] + torque*(Ir/tang) - pend_fric*y[1];
      else if (comp == 2) return 0.0;
      else if (comp == 3) return 0.0;
      else{
        std::cerr << "No more dofs." << '\n';
        return 0.0;
      }
    }
    else if(phi_max < d && d <= phi/2.0){
      tag_phase = 2;                            // First uncoupled phase (B1)
      if (comp == 0) return y[1];
      else if (comp == 1) return -omega*omega*y[0] - pend_fric*y[1];
      else if (comp == 2) return y[3];
      else if (comp == 3) return -pend_fric*y[3]+torque;
      else{
        std::cerr << "No more dofs." << '\n';
        return 0.0;
      }
    }
    else if(phi/2.0 < d && d <= (phi_max + phi/2.0)){
      tag_phase = 3;                            // Second coupled phase (A2)
      double tang = tan(alpha);
      double tang2 = tang*tang;
      double I = (pend_inertia*tang2 + wheel_inertia)/tang2;
      double Tr = torque*(wheel_inertia/(I*tang));
      double Io = pend_inertia/I;
      double Ir = wheel_inertia/I;
      if (comp == 0) return y[1];
      else if (comp == 1) return -Io*omega*omega*y[0] - torque*(Ir/tang) - pend_fric*y[1];
      else if (comp == 2) return 0.0;
      else if (comp == 3) return 0.0;
      else{
        std::cerr << "No more dofs." << '\n';
        return 0.0;
      }
    }
    else if((phi_max + phi/2.0) < d && d <= phi){
      tag_phase = 4;                            // Second uncoupled phase (B2)
      if (comp == 0) return y[1];
      else if (comp == 1) return -omega*omega*y[0] - pend_fric*y[1];
      else if (comp == 2) return y[3];
      else if (comp == 3) return -pend_fric*y[3]+torque;
      else{
        std::cerr << "No more dofs." << '\n';
        return 0.0;
      }
    }
    else{
      std::cerr << "Outside allowed range of Cogwheel." << '\n';
      return 0.0;
    }
  }
  else{
    std::cerr << "Outside allowed range of Pendullum." << '\n';
    return 0.0;
  }
}
/*******************************************************************************
                  INITIALIZING MEMBER FUNCTIONS: DECLARATION
*******************************************************************************/
void Clock::init_design_param(double phi, double r, double a, double teeth, double fric){
  phi_max = phi;
  alpha = a;
  rho = r;
  num_teeth = teeth;
  int_fric = fric;
}
/******************************************************************************/
void Clock::init_pend_mec_param(double w, double I, double fric){
  omega = w;
  pend_inertia = I;
  pend_fric = fric;
}
/******************************************************************************/
void Clock::init_wheel_mec_param(double tau, double I, double fric){
  torque = tau;
  wheel_inertia = I;
  wheel_fric = fric;
}
/******************************************************************************/
