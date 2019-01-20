#include "Punchet1.h"

int main(void){
  Punchet clock1;
  clock1.init_punchet();
  double Phi = M_PI/clock1.num_teeth;

  clock1.now_init_cond[0] = 0.9*clock1.theta_max;
  clock1.now_init_cond[1] = (M_PI/180.0)*0.0*clock1.theta_max*clock1.omega;
  clock1.now_init_cond[2] = 0.0;
  clock1.now_init_cond[3] = 0.0;
  clock1.now_init_cond[4] = 0.0;

  /*
  for(int ii = 1; ii <= 2; ii++){
    printf("FINDING TIME OF COUPLING\n");
    clock1.t_coupling(clock1.now_phase,clock1.coup_num);
    clock1.mapping();
    printf("----------------------------------------------------------------\n");
    clock1.print_state();
    clock1.print_fut_state();
    clock1.update_init_cond();
    printf("****************************************************************\n");
  }
  */


  for(int ii = 1; ii<= 400; ii++){
    clock1.mapping();
    clock1.print_motion(1000);
    clock1.update_init_cond();
  }


  /*
  Print_title();
  for(int ii = 1; ii<= 300; ii++){
    clock1.print_state();
    clock1.mapping();
    //printf("-----------------------------------------------------------------------------------------\n");
    clock1.update_init_cond();
  }
  */

  /*
  printf("%4.7f\n",clock1.theta_max);
  for(int ii = 1; ii <= 4; ii++){
    double tc = clock1.t_coupling(clock1.now_phase,ii);
    double fc = clock1.coupling_func(tc,clock1.now_phase,ii);
    double theta = clock1.theta(tc,clock1.now_phase);
    printf("%d\t%d\t%4.7f\t%4.7e\t%4.7f\t",clock1.now_phase,ii,tc,fc,theta);
    printf("%d\n",clock1.coup_index(tc,ii));
  }
  double t = clock1.t_transit();
  printf("Coup_num = %d\tTime coup. = %4.7f\n",clock1.coup_num,t);
  */

  /*
  printf("%4.7f\n",clock1.theta_max);
  for(int ii = 1; ii<= 3; ii++){
    clock1.print_state();
    double tr = clock1.t_coupling(clock1.now_phase,clock1.coup_num);
    double tr1 = clock1.t_coupling(clock1.now_phase,clock1.coup_num+1);
    double fr = clock1.coupling_func(tr,clock1.now_phase,clock1.coup_num);
    double fr1 = clock1.coupling_func(tr1,clock1.now_phase,clock1.coup_num+1);
    double theta = clock1.theta(tr,clock1.now_phase);
    double theta1 = clock1.theta(tr1,clock1.now_phase);
    printf("-----------------------------------------------------------------\n");
    printf("%4.7f\t%4.7f\t%4.7f\n",tr,fr,theta);
    printf("%4.7f\t%4.7f\t%4.7f\n",tr1,fr1,theta1);
    printf("-----------------------------------------------------------------\n");
    clock1.mapping();
    clock1.update_init_cond();
    clock1.print_fut_state();
    printf("*****************************************************************\n");

  }
  */

  /*
  for(int ii = 0; ii <= 1000; ii++){
    double t = clock1.now_init_cond[4] + 0.001*ii;
    double theta = clock1.theta(t,2);
    double theta_dot = clock1.theta_dot(t,2);
    double phi = clock1.phi(t,2);
    double phi_dot = clock1.phi_dot(t,2);
    printf("%4.9f\t%4.9f\t%4.9f\t%4.9f\t%4.9f\n",t,phi,phi_dot,theta,theta_dot);
  }

  clock1.mapping();
  clock1.update_init_cond();

  for(int ii = 0; ii <= 100; ii++){
    double t = clock1.now_init_cond[4] + 0.001*ii;
    double theta = clock1.theta(t,clock1.now_phase);
    double phi = clock1.phi(t,clock1.now_phase);
    printf("%4.9f\t%4.9f\t%4.9f\n",t,phi,theta);
  }
  */

  /*
  clock1.mapping();

  printf("CURRENT STATE\n");
  clock1.print_state();
  printf("FUTURE STATE\n");
  clock1.print_fut_state();
  */

  /*
  clock1.coup_num = 1;
  int num = clock1.coup_num;



  double t0 = clock1.now_init_cond[4] = 0.0;
  double t_max = clock1.t_coupling(1,num);
  double t_max1 = clock1.t_transit();

  //std::cout << t_max << "\t" << t_max1 << "\n";



  double t = t0;
  for(int ii = 0; t <= t_max; ii++){
    t += 0.01;

    double theta = clock1.theta(t,1);
    double phi = clock1.phi(t,1);

    std::printf("%4.9f\t%4.9f\t%4.9f\n",t,phi,theta);
  }
  */

  clock1.kill_punchet();
  return 0;
}
