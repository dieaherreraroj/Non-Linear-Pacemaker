#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <vector>

struct Punchet{
/*******************************************************************************
                        DESING PARAMETERS OF ESCAPEMENT
*******************************************************************************/
  double alfa;                             // Constraint parameter
  double rho;                              // radii ratio
  double num_teeth;                        // # of teeth of escapement wheel
  double theta_max;                        // Top angle of Penulum (coup. mot.)
  double torque;                           // External Torque for esc. wheel
/*******************************************************************************
                      MECHANICAL PARAMETERS OF ESCAPEMENT
*******************************************************************************/
  double omega;                            // Natural Frequency of Pendulum
  double Ia;                               // Inertia Moment of Pendulum
  double Iw;                               // Inertia Moment of escapement wheel
  double Itot;                             // Coupled Moment of Inertia
  double gamma;                            // Friction Coefficient of Pendulum
  double gamma_w;                          // Friction Coefficient of esc. wheel
  double Gamma_tot;                        // Friction coeff. coupled motion
  int coup_num;                            // Type of coupling indicator
/*******************************************************************************
                 INITIAL-FINAL CONDITIONS FOR UPCOMING MOTION
*******************************************************************************/
  double *now_init_cond;                   // Contains init. cond. and init time
  int now_phase;                           // Indicates current phase of motion
  double *next_init_cond;                  // Contains final cond. and init time
  int next_phase;                          // Indicates next phase of motion
  int next_coup_num;                       // Type of coupling indicator
/*******************************************************************************
                             DEFAULT INTIALIZER
*******************************************************************************/
  Punchet(){
    alfa = tan(30.0*M_PI/180.0);            // = 0.0175
    rho = 1.2;
    num_teeth = 30.0;
    theta_max = M_PI/1080;   // < PHI/(4.0*alfa)
    omega = 2*M_PI;
    gamma = 0.2;
    gamma_w = 0.9*gamma;
    Gamma_tot = gamma + alfa*alfa*gamma_w/rho;
    Ia = 1.0;
    Iw = 0.9*Ia;
    Itot = Ia + Iw*alfa*alfa/rho;
    torque = (1.7)*omega*omega*theta_max*Ia/alfa;
    now_phase = 1;
    next_phase = 2;
    coup_num = 1;
  }
/*******************************************************************************
                 FUNCTIONS FOR DIFFERENT DEGREES OF FREEDOM
*******************************************************************************/
  double theta(double t, int tag);         // Angle of anchor and pendulum
  double theta_dot(double t, int tag);     // Velocity of anchor and pendulum
  double phi(double t, int tag);           // Angle of escapement wheel
  double phi_dot(double t, int tag);       // Velocity of escapement wheel
/*******************************************************************************
                           AUXILIAR FUNCTIONS
*******************************************************************************/
  void init_punchet();
  void kill_punchet();
  void update_init_cond();
  void print_state();
  void print_fut_state();
  void print_motion(int num);
  double coupling_func(double t, int tag, int num);
  double t_coupling(int tag, int num);
  int coup_index(double t, int num);
  double t_transit();
  void mapping();
};
/*******************************************************************************
                            ALLOCATION OF MEMORY
*******************************************************************************/
void Punchet::init_punchet(){
  now_init_cond = (double*) calloc(5,sizeof(double));
  next_init_cond = (double*) calloc(5,sizeof(double));
}
void Punchet::kill_punchet(){
  free(now_init_cond);
  free(next_init_cond);
}
void Punchet::update_init_cond(){
  now_phase = next_phase;
  coup_num = next_coup_num;
  for(int ii = 0; ii<5; ii++)
    now_init_cond[ii] =  next_init_cond[ii];
}
/*******************************************************************************
                             PRINTING FUNCTIONS
*******************************************************************************/
void Punchet::print_state(){
  double psi = M_PI/num_teeth;
  double Theta = now_init_cond[0]/theta_max;
  double Theta_dot = now_init_cond[1]/(omega*theta_max);
  double Phi = now_init_cond[2]/psi;
  double Phi_dot = now_init_cond[3]/(gamma_w*psi);
  double t = now_init_cond[4];
  printf("%d\t",now_phase);
  printf("%4.9f\t%4.9f\t%4.9f\t%4.9f\t%4.9f\t",t,Phi,Phi_dot,Theta,Theta_dot);
  /*
  for(int ii = 4; ii >= 0; ii--)
    printf("%4.10f\t",now_init_cond[ii]);
  */
  printf("%d\n",coup_num);

}
void Punchet::print_fut_state(){
  printf("%d\t",next_phase);
  for(int ii = 4; ii >= 0; ii--)
    printf("%4.10f\t",next_init_cond[ii]);
  printf("%d\n",coup_num);
}
void Punchet::print_motion(int num){
  double ti = now_init_cond[4];
  double tf = next_init_cond[4];
  double dt = (tf-ti)/num;
  double Theta = 0.0;
  double Theta_dot = 0.0;
  double Phi = 0.0;
  double Phi_dot = 0.0;
  double t = ti;
  double psi = M_PI/num_teeth;
  for(int ii = 0; t <= tf; ii++){
    t += dt;
    double Theta = theta(t,now_phase)/theta_max;
    double Theta_dot = theta_dot(t,now_phase)/(omega*theta_max);
    double Phi = phi(t,now_phase)/psi;
    double Phi_dot = phi_dot(t,now_phase)/(gamma_w*psi);
    printf("%4.9f\t%4.9f\t%4.9f\t%4.9f\t%4.9f\n",t,Phi,Phi_dot,Theta,Theta_dot);
  }
}
/*******************************************************************************
                    ANALYTICAL SOLUTIONS OF EOM FOR ANCHOR
*******************************************************************************/
double Punchet::theta(double t, int tag){
  if(tag == 1){
    // Parameters derived from analytical integration
    double b = gamma/Ia;
    double w2 = omega*omega-(b/2.0)*(b/2.0);
    double w = sqrt(fabs(w2));
    double B = now_init_cond[0];
    double A = (now_init_cond[1]+0.5*b*B)/w;
    double tp = t - now_init_cond[4];
    // Case of damped oscillatory motion
    if(w2 > 0.0) return exp(-0.5*b*tp)*(A*sin(w*tp)+B*cos(w*tp));
    // Case of overdamped motion
    else if(w2 < 0.0) return exp(-0.5*b*tp)*(A*sinh(w*tp)+B*cosh(w*tp));
    // Case of critically damped motion
    else{
      A = now_init_cond[1] + 0.5*b*B;
      return (B+A*tp)*exp(-0.5*b*tp);
    }
  }
  else if (tag == 2){
    // Parameters derived from analytical integration
    double b = gamma/Itot;
    double k = Ia/Itot;
    double tau = pow(-1,coup_num)*alfa*torque/rho;
    double d0 = now_init_cond[0] - tau/(Ia*omega*omega);
    double w2 = k*omega*omega-(b/2.0)*(b/2.0);
    double w = sqrt(fabs(w2));
    double B = d0;
    double A = (now_init_cond[1]+0.5*b*B)/w;
    double tp = t - now_init_cond[4];
    // Change of variable according to external torque
    double delta = 0;
    if(w2 > 0.0) delta = exp(-0.5*b*tp)*(A*sin(w*tp)+B*cos(w*tp));
    // Case of overdamped motion
    else if(w2 < 0.0) delta = exp(-0.5*b*tp)*(A*sinh(w*tp)+B*cosh(w*tp));
    // Case of critically damped motion
    else{
      A = now_init_cond[1] + 0.5*b*B;
      delta = (B+A*tp)*exp(-0.5*tp);
    }
    return delta + tau/(Ia*omega*omega);
  }
  else{
    std::cerr << "There is not another phase of motion" << '\n';
    return -1000;
  }
}
double Punchet::theta_dot(double t, int tag){
  if(tag == 1){
    // Parameters derived from analytical integration
    double b = gamma/Ia;
    double w2 = omega*omega-(b/2.0)*(b/2.0);
    double w = sqrt(fabs(w2));
    double B = now_init_cond[0];
    double A = (now_init_cond[1]+0.5*b*B)/w;
    double tp = t - now_init_cond[4];
    // Case of damped oscillatory motion
    if(w2 > 0.0) return exp(-0.5*b*tp)*((w*A-0.5*b*B)*cos(w*tp)-(w*B+0.5*b*A)*sin(w*tp));
    // Case of overdamped motion
    else if(w2 < 0.0) return exp(-0.5*b*tp)*((w*A-0.5*b*B)*cosh(w*tp)+(w*B-0.5*b*A)*sinh(w*tp));
    // Case of critically damped motion
    else{
      A = now_init_cond[1] + 0.5*b*B;
      return (A-0.5*b*B-0.5*b*A*tp)*exp(-0.5*b*tp);
    }
  }
  else if(tag == 2){
    // Parameters derived from analytical integration
    double b = gamma/Itot;
    double k = Ia/Itot;
    double tau = pow(-1,coup_num)*alfa*torque/rho;
    double d0 = now_init_cond[0] - tau/(omega*omega*alfa);
    double w2 = k*omega*omega-(b/2.0)*(b/2.0);
    double w = sqrt(fabs(w2));
    double B = d0;
    double A = (now_init_cond[1]+0.5*b*B)/w;
    double tp = t - now_init_cond[4];
    // Change of variable according to external torque
    double delta = 0;
    if(w2 > 0.0) delta = exp(-0.5*b*tp)*((w*A-0.5*b*B)*cos(w*tp)-(w*B+0.5*b*A)*sin(w*tp));
    // Case of overdamped motion
    else if(w2 < 0.0) delta = exp(-0.5*b*tp)*((w*A-0.5*b*B)*cosh(w*tp)+(w*B-0.5*b*A)*sinh(w*tp));
    // Case of critically damped motion
    else{
      A = now_init_cond[1] + 0.5*b*B;
      delta = (A-0.5*b*B-0.5*b*A*tp)*exp(-0.5*b*tp);
    }
    return delta;
  }
  else{
    std::cerr << "There is not another phase of motion" << '\n';
    return -1000;
  }
}
/*******************************************************************************
                    ANALYTICAL SOLUTIONS OF EOM FOR ANCHOR
*******************************************************************************/
double Punchet::phi(double t, int tag){
  if(tag == 1){
    // Parameters derived from analytical integration
    double b = gamma_w/Iw;
    double f = torque/Iw;
    double tp = t - now_init_cond[4];
    return (f/b)*tp+now_init_cond[2]-((f/b-now_init_cond[3])*(1.0-exp(-b*tp)))/b;
  }
  else if (tag == 2) {
    return pow(-1,coup_num)*(alfa/rho)*theta(t,tag)+(coup_num*0.5)*M_PI/num_teeth;
  }
  else{
    std::cerr << "There is not another phase of motion" << '\n';
    return -1000;
  }
}
double Punchet::phi_dot(double t, int tag){
  if(tag == 1){
    // Parameters derived from analytical integration
    double b = gamma_w/Iw;
    double f = torque/Iw;
    double tp = t - now_init_cond[4];
    return (f/b)-((f/b-now_init_cond[3])*(exp(-b*tp)));
  }
  else if (tag == 2){
    return pow(-1,coup_num)*(alfa/rho)*theta_dot(t,tag);
  }
  else{
    std::cerr << "There is not another phase of motion" << '\n';
    return -1000;
  }}
/*******************************************************************************
                      COUPLING - DECOUPLING FUNCTIONS
*******************************************************************************/
double Punchet::coupling_func(double t, int tag, int num){
  if(tag == 1){
    double Phi = M_PI/num_teeth;                      // = 0.1
    return pow(-1,num)*(rho/alfa)*(phi(t,tag)-0.5*Phi*num)-theta(t,tag);
  }
  else if (tag == 2){
    return theta(t,tag)-pow(-1,coup_num)*theta_max;
  }
  else{
    std::cerr << "There is not another phase of motion" << '\n';
    return -1000;
  }
}
double Punchet::t_coupling(int tag, int num){
  // Implementation of Regula Falsi algorithm.
  double t1 = now_init_cond[4];             // First guess root
  double t2 = t1;                           // Second guess root
  double tr = t1;                           // Time root
  // Determine existance of root by continuity
  double f1 = coupling_func(t1,tag,num);
  double f2 = coupling_func(t2,tag,num);
  for(int ii = 0; ii < 100000; ii++){
    if(f1*f2 < 0.0)
      break;
    else{
      t2 += 0.01;
      f2 = coupling_func(t2,tag,num);
    }
  }
  if(t1 < t2){
    double fr = -1000;
    for(int ii = 1; ii <= 10000; ii++){
      f1 = coupling_func(t1,tag,num);
      f2 = coupling_func(t2,tag,num);
      // Regula Falsi Formula
      tr = (f2*t1-f1*t2)/(f2-f1);
      // Check whether root or not
      fr = coupling_func(tr,tag,num);
      if(std::fabs(fr) < 1e-13)
        break;
      else{
        if(fr*f1 > 0.0)
          t1 = tr;
        else
          t2 = tr;
      }
    }
    if(std::fabs(fr) < 1e-13)
      return tr;
    else{
      std::cerr << "False root found, try again, idiot." << '\n';
      return -10000;
    }
  }
  else{
    std::cerr << "No root found, try again" << '\n';
    return -10000;
  }
}
int Punchet::coup_index(double t, int num){
  double c = pow(-1,num+1)*(theta(t,now_phase)-pow(-1,num)*theta_max);
  if(c > 0.0)
    return 1;
  else
    return 0;
}
double Punchet::t_transit(){
  // Variable for storing time of coupling
  double tc = -1000;
  // Carry out mapping process (for initial conditions)
  if(now_phase == 1){                        // From decoupling to coupling
    for(int ii = 0; ii <= 10; ii++){
      int tag = coup_num+ii;
      tc = t_coupling(now_phase,tag);
      int index = coup_index(tc,tag);
      if(index == 1){
        coup_num = tag;
        break;
      }
      else
        tc = -1000;
    }
    return tc;
    /*
    // Check wich type of coupling is next
    double t1 = t_coupling(now_phase,coup_num);
    double t2 = t_coupling(now_phase,coup_num+1);
    // Compute the possible anchor angles of coupling
    double theta1 = theta(t1,now_phase);
    double theta2 = theta(t2,now_phase);
    // Check correct line of coupling (see document)
    if((t1 < t2) && (pow(-1,coup_num+1)*(theta1 - pow(-1,coup_num)*theta_max) > 0.0)){
      tc = t1;
      attempts = 1;
      //printf("hello sweety.\n");
    }
    else if (pow(-1,coup_num+2)*(theta2 - pow(-1,coup_num+1)*theta_max) > 0.0){
      tc = t2;
      //printf("Hello darling.\n");
      coup_num = coup_num + 1;
      attempts = 1;
    }
*/
  }
  else if(now_phase == 2)                   // From coupling to decoupling
    tc = t_coupling(now_phase,coup_num);
  else{
    std::cerr << "There is not another phase of motion" << '\n';
    tc = -2000;
  }
  return tc;
}
void Punchet::mapping(){
  //printf("%d\n",now_phase);
  double tc = t_transit();
  //printf("%4.7f\n",tc);
  if(tc != -1000 && tc > now_init_cond[4]){
    // From decoupling to coupling
    if(now_phase == 1){
      next_phase = 2;
      // Store future initial conditions
      next_init_cond[0] = theta(tc,now_phase);
      next_init_cond[1] = theta_dot(tc,now_phase);
      next_init_cond[2] = phi(tc,now_phase);
      // Careful with impact boundary conditions (!!!)
      next_init_cond[3] = pow(-1,coup_num)*(alfa/rho)*theta_dot(tc,now_phase);
      next_init_cond[4] = tc;
      next_coup_num = coup_num;
    }
    // From Coupling to decoupling
    else if (now_phase == 2) {
      next_phase = 1;
      // Store future initial conditions
      next_init_cond[0] = theta(tc,now_phase);
      next_init_cond[1] = theta_dot(tc,now_phase);
      next_init_cond[2] = phi(tc,now_phase);
      next_init_cond[3] = phi_dot(tc,now_phase);
      next_init_cond[4] = tc;
      next_coup_num = coup_num + 1;                                 // To avoid repetition
    }
    else
      std::cerr << "There is not another phase of motion." << '\n';
  }
  else if(tc == -1000){
    std::cerr << "Coupling number top in t_transit too small, idiot." << '\n';
  }
  else if(tc == -2000)
    std::cerr << "There is not another phase of motion." << "\n";
  else
    std::cerr << "Check errors in finding root time, idiot." << '\n';
}

void Print_title(){
  printf("T\tTime_transit\tPhi (wheel)\tPhi_d(wheel)\tTheta (anc) \tTheta_d(anc)\tC\n");
}
