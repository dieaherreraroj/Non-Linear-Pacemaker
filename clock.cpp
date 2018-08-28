#include "Punchet.h"

int main(void){

  // Declaration of variables
  Punchet grandpa;
  std::vector<double> data;
  double PHI = (2.0*M_PI/grandpa.num_teeth);
  // Request memory
  grandpa.init_punchet();

  // Set intitial conditions for motion;
  //grandpa.num_cycle = 6;
  grandpa.now_init_cond[0] = 0.0;
  grandpa.now_init_cond[1] = 0.2;
  grandpa.now_init_cond[2] = 0.0;
  grandpa.now_init_cond[3] = 0.0;
  grandpa.now_init_cond[4] = 0.0;

  //grandpa.print_table(200.0);
  grandpa.num_solve(200.0,data);

  // Return memory
  grandpa.kill_punchet();
  return 0;
}
