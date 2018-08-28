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
  double num_teeth;                        // # of teeth of escapement wheel
  double gamma;                            // Friction Coefficient of coup. mot.
  double torque;                           // External Torque for esc. wheel
  double theta_max;                        // Top angle of Penulum (coup. mot.)
/*******************************************************************************
                      MECHANICAL PARAMETERS OF ESCAPEMENT
*******************************************************************************/
  double It;                               // Inertia Moment of Pendulum
  double Ip;                               // Inertia Moment of escapement wheel
  double Itot;                             // Coupled Moment of Inertia
  double gammat;                           // Friction Coefficient of Pendulum
  double gammap;                           // Friction Coefficient of esc. wheel
  double omega;                            // Natural Frequency of Pendulum
  int num_cycle;                        // Number of tooth in current cycle
/*******************************************************************************
                             DEFAULT INTIALIZER
*******************************************************************************/
  Punchet(){
/******************************************************************************/
    alfa = tan(5.0*M_PI/180.0);
    num_teeth = 40.0;
    gamma = 0.2;
    theta_max = 0.01*(2.0*M_PI/num_teeth)/alfa;   // < PHI/(4.0*alfa)
    omega = 2.0*M_PI;
/******************************************************************************/
    It = 1.0;
    Ip = 0.1*It;
    Itot = It + Ip/(alfa*alfa);
    gammat = 0.85*gamma;
    gammap = 0.01*gammat;
    torque = 1.10*omega*omega*theta_max*It*alfa;
    num_cycle = 0;
  }
/*******************************************************************************
                 INITIAL-FINAL CONDITIONS FOR UPCOMING MOTION
*******************************************************************************/
  double *now_init_cond;                   // Contains init. cond. and init time
  int now_phase;                           // Indicates current phase of motion
  double *next_init_cond;                  // Contains final cond. and init time
  int next_phase;                          // Indicates next phase of motion
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
  void check_now_phase();
  void check_num_cycle();
  double change_phase(double t, int tag1, int tag2);
  double time_change_phase(int tag1, int tag2, double t_top);
  void init_next_phase(double t_top);      // Here is the real magic
  void print_phase(std::vector<double> & data);
  void update_init_cond();
  void num_solve(double t_top, std::vector<double> & data);
  void print_table(double t_top);
};
/*******************************************************************************
           FUNCTIONS FOR DIFFERENT DEGREES OF FREEDOM: IMPLEMENTATIONS
*******************************************************************************/
double Punchet::theta(double t, int tag){
  if(tag == 3){
    // Parameters derived from analytical integration
    double b = gammat/It;
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
  else if((tag == 1) || (tag == 2)){
    // Parameters derived from analytical integration
    double b = gamma/Itot;
    double k = It/Itot;
    double d0 = now_init_cond[0] + pow(-1.0,tag)*torque/(k*omega*omega*alfa);
    double w2 = k*omega*omega-(b/2.0)*(b/2.0);
    double w = sqrt(fabs(w2));
    double B = d0;
    double A = (now_init_cond[1]+0.5*b*B)/sqrt(fabs(w2));
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
    return delta - pow(-1.0,tag)*torque/(k*omega*omega*alfa);
  }
  else{
    std::cerr << "There is not another phase of motion" << '\n';
    return 0.0;
  }
}
/******************************************************************************/
double Punchet::theta_dot(double t, int tag){
  if(tag == 3){
    // Parameters derived from analytical integration
    double b = gammat/It;
    double w2 = omega*omega-(b/2.0)*(b/2.0);
    double w = sqrt(fabs(w2));
    double B = now_init_cond[0];
    double A = (now_init_cond[1]+0.5*b*B)/w;
    double tp = t - now_init_cond[4];
    // Case of damped oscillatory motion
    if(w2 > 0.0) return exp(-0.5*b*tp)*((w*A-0.5*b*B)*cos(w*tp)-(w*B+0.5*A)*sin(w*tp));
    // Case of overdamped motion
    else if(w2 < 0.0) return exp(-0.5*b*tp)*((w*A-0.5*b*B)*cosh(w*tp)+(w*B-0.5*A)*sinh(w*tp));
    // Case of critically damped motion
    else{
      A = now_init_cond[1] + 0.5*b*B;
      return (A-0.5*b*B-0.5*b*A*tp)*exp(-0.5*b*tp);
    }
  }
  else if(tag == 1 || tag == 2){
    // Parameters derived from analytical integration
    double b = gamma/Itot;
    double k = It/Itot;
    double d0 = now_init_cond[0] + pow(-1.0,tag)*torque/(k*omega*omega*alfa);
    double w2 = (k*omega)*(k*omega)-(b/2.0)*(b/2.0);
    double w = sqrt(fabs(w2));
    double B = d0;
    double A = (now_init_cond[1]+0.5*b*B)/sqrt(fabs(w2));
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
    return 0.0;
  }
}
/******************************************************************************/
double Punchet::phi(double t, int tag){
  if(tag == 3){
    // Parameters derived from analytical integration
    double b = gammap/Ip;
    double f = torque/Ip;
    double tp = t - now_init_cond[4];
    return (f/b)*tp+now_init_cond[2]-((f/b-now_init_cond[3])*(1.0-exp(-b*tp)))/b;
  }
  else if(tag == 1){
    double PHI = 2.0*M_PI/num_teeth;
    return Punchet::theta(t,tag)*alfa+PHI*num_cycle;
  }
  else if(tag == 2){
    double PHI = 2.0*M_PI/num_teeth;
    return -Punchet::theta(t,tag)*alfa+PHI*num_cycle+0.5*PHI;
  }
  else{
    std::cerr << "There is not another phase of motion" << '\n';
    return 0.0;
  }
}
/******************************************************************************/
double Punchet::phi_dot(double t, int tag){
  if(tag == 3){
    // Parameters derived from analytical integration
    double b = gammap/Ip;
    double f = torque/Ip;
    double tp = t - now_init_cond[4];
    return (f/b)-((f/b-now_init_cond[3])*(exp(-b*tp)));
  }
  else if(tag == 1){
    double PHI = 2.0*M_PI/num_teeth;
    return Punchet::theta_dot(t,tag)*alfa;
  }
  else if(tag == 2){
    double PHI = 2.0*M_PI/num_teeth;
    return -Punchet::theta_dot(t,tag)*alfa;
  }
  else{
    std::cerr << "There is not another phase of motion" << '\n';
    return 0.0;
  }
}
/*******************************************************************************
                           AUXILIAR FUNCTIONS
*******************************************************************************/
void Punchet::init_punchet(){
  now_init_cond = (double*) calloc(5,sizeof(double));
  next_init_cond = (double*) calloc(5,sizeof(double));
}
/******************************************************************************/
void Punchet::kill_punchet(){
  free(now_init_cond);
  free(next_init_cond);
}
/******************************************************************************/
void Punchet::check_now_phase(){
  double PHI = 2.0*M_PI/num_teeth;
  double dtic = 0.0, tic = 0.0, dtac = 0.0, tac = 0.0;
  dtic = (now_init_cond[2]-now_init_cond[0]*alfa)/PHI;
  tic = dtic - floor(dtic);
  dtac = (now_init_cond[2]+now_init_cond[0]*alfa)/PHI-0.5;
  tac = dtac - floor(dtac);
  if(fabs(tic) < 1e-10){
    if(now_init_cond[0] <= theta_max)
      now_phase = 1;
    else
      now_phase = 3;
  }
  else{
    if(fabs(tac) < 1e-10){
      if(-theta_max <= now_init_cond[0])
        now_phase = 2;
      else
        now_phase = 3;
    }
    else
      now_phase = 3;
  }
}
/******************************************************************************/
void Punchet::check_num_cycle(){
  double PHI = 2.0*M_PI/num_teeth;
  if(PHI*(0.5+num_cycle) <= next_init_cond[2]){
    //printf("next tooth   ");
    num_cycle++;
  }
  else{
    //printf("stay tooth   ");
  }
}
/******************************************************************************/
double Punchet::change_phase(double t, int tag1, int tag2){
  double PHI = 2.0*M_PI/num_teeth;
  if((tag1 == 1) && (tag2 == 3)) return theta(t,tag1)-theta_max;
    else if((tag1 == 2) && (tag2 == 3)) return theta(t,tag1)+theta_max;
      else if((tag1 == 3) && (tag2 == 1)){
        return theta(t,tag1)*alfa-phi(t,tag1)+num_cycle*PHI;
      }
      else if((tag1 == 3) && (tag2 == 2)){
        return theta(t,tag1)*alfa+phi(t,tag1)-(num_cycle+0.5)*PHI;
      }
  else{
    //std::cerr << "Trying to keep same phase or unidentified phase" << '\n';
    return 0.0;
  }
}
/******************************************************************************/
double Punchet::time_change_phase(int tag1, int tag2, double t_top){
  double ti = now_init_cond[4], tf, tr;
  double fi = change_phase(ti,tag1,tag2), ff = fi, fr;
  // Loof for next seed for Regula Falsi
  for(int ii = 0; ii <= 10000; ii++){
    tf = ti + ii*0.0001;
    ff = change_phase(tf,tag1,tag2);
    if(ff*fi > 0.0)
      continue;
    else{
      if(ff*fi < 0.0)
        break;
      else{
        if(ff == 0.0)
          tr = tf;
        else
          tr = ti;
        break;
      }
    }
  }
  int ii;
  // Apply Regula Falsi Algorythm
  for(ii = 1; ii <= 10000; ii++){
    fi = change_phase(ti,tag1,tag2);
    ff = change_phase(tf,tag1,tag2);
    tr = tf - (ff*(ti-tf))/(fi-ff);
    fr = change_phase(tr,tag1,tag2);
    if(fabs(fr) < 1e-14){
      break;
    }
    else{
      if(fr*fi > 0.0)
        ti = tr;
      else
        tf = tr;
    }
  }
  // Check whether tr is an actual root
  fr = change_phase(tr,tag1,tag2);
  if(ii == 10000){
    if(fabs(fr) > 1e-14){
      std::cerr << "Need more steps or there is no root " << '\n';
      return t_top;
    }
    else{
      std::cerr << "Finished iterations " << '\n';
      return tr;
    }
  }
  else
    //std::cout << "Hello Diego!!!" << tr <<'\n';
    return tr;
}
/******************************************************************************/
void Punchet::init_next_phase(double t_top){
  if((now_phase == 1) || (now_phase == 2)){
    // Looking for uncoupled motion
    double tr = time_change_phase(now_phase,3,t_top);
    // Idenify actual end of fase (including ending of motion record)
    double t = fmin(tr,t_top);
    // Store next initial conditions of clock motion
    next_init_cond[4] = t;
    next_init_cond[0] = theta(t,now_phase);
    next_init_cond[1] = theta_dot(t,now_phase);
    next_init_cond[2] = phi(t,now_phase);
    next_init_cond[3] = phi_dot(t,now_phase);
    // Determine next phase of motion
    if(t == tr)
      next_phase = 3;
    else
      next_phase = now_phase;
  }
  else if(now_phase == 3){
    double t1r = time_change_phase(now_phase,1,t_top);
    double t2r = time_change_phase(now_phase,2,t_top);
    double tp = 0.0;
    if((theta(t1r,now_phase) <= theta_max) && (-theta_max  <= theta(t2r,now_phase)))
      tp = fmin(t1r,t2r);
    else{
      if(theta(t1r,now_phase) <= theta_max) tp = t1r;
        else if(-theta_max  <= theta(t2r,now_phase)) tp = t2r;
      else
        tp = t_top;
    }
    double t = fmin(tp,t_top);
    // Store next initial conditions of clock motion
    next_init_cond[4] = t;
    next_init_cond[0] = theta(t,now_phase);
    next_init_cond[1] = theta_dot(t,now_phase);
    next_init_cond[2] = phi(t,now_phase);
    // Apply boundary conditions of collision
    // And determine next phase
    if(t == tp){
      next_init_cond[3] = alfa*next_init_cond[1];
      if(tp == t1r) next_phase = 1;
        else if(tp == t2r) next_phase = 2;
      else
        next_phase = now_phase;
    }
    else{
      next_init_cond[3] = phi_dot(t,now_phase);
      next_phase = now_phase;
    }
  }
  else{
    std::cerr << "Bad Phase of motion, check now_phase" << '\n';
  }
}
/******************************************************************************/
void Punchet::print_phase(std::vector<double> & data){
  double PHI = 2.0*M_PI/30.0;
  if((0<now_phase) && (now_phase<4)){
    if((0<next_phase) && (next_phase<4)){
      double t = now_init_cond[4];
      double dt = 1e-4;
      for(int ii = 0; t <= next_init_cond[4]; ii++){
        t = now_init_cond[4] + ii*dt;
        double y = theta(t,now_phase);
        data.push_back(y);
        double vy = theta_dot(t,now_phase);
        double x = phi(t,now_phase);
        double vx = phi_dot(t,now_phase);
        y *= 180.0/M_PI;
        x *= 1.0/PHI;
        printf("%1d\t %4.7e\t %4.7e\t %4.7e\t %4.7e\t %4.7e\n",now_phase,t,y,vy,x,vx);
      }
    }
    else{
      std::cerr << "Not updated next initial conditions" << '\n';
    }
  }
  else{
    std::cerr << "Not updated now initial conditions" << '\n';
  }
}
/******************************************************************************/
void Punchet::update_init_cond(){
  for(int ii = 0; ii<5; ii++)
    now_init_cond[ii] = next_init_cond[ii];
  now_phase = next_phase;
}
/******************************************************************************/
void Punchet::num_solve(double t_top, std::vector<double> & data){
  // Check initial phase of motion
  check_now_phase();
  check_num_cycle();
  // Compute phase transtions in some interval [0.0:t_top]
  while(next_init_cond[4] < t_top){
    // Set up next phase of motion
    init_next_phase(t_top+1.0);
    // Print data and record it (timestep == 1e-4)
    if(now_phase == 1)
    print_phase(data);
    // Update initial conditions: next -> now
    update_init_cond();
    // Update phase of motion
    check_num_cycle();
  }
}
/******************************************************************************/
void Punchet::print_table(double t_top){
  double PHI = 2.0*M_PI/num_teeth;
  check_now_phase();
  check_num_cycle();
  printf("NP\tNxP\t#T\t      tr\t   vi\t   vf\t\tCF\n");
  while(next_init_cond[4] < t_top){
    // Set up next phase of motion
    init_next_phase(t_top+1.0);
    printf("%1d\t ",now_phase);
    printf("%1d\t ",next_phase);
    printf("%1d\t ",num_cycle);
    double tr = next_init_cond[4];
    double phi_init = now_init_cond[2];
    double phi_next = next_init_cond[2];
    printf("%4.7e\t %4.3f\t %4.3f\t ",tr,phi_init/PHI,phi_next/PHI);
    printf("(%4.7e)\n",change_phase(tr,now_phase,next_phase));
    update_init_cond();
    check_num_cycle();
  }
}
