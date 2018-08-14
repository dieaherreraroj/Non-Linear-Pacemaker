#include "Punchet.h"

int main(void){
  // Declaration of variables
  Punchet grandpa; std::vector<double> data;
  // Request memory
  grandpa.init_punchet();
  // Set intitial conditions for motion;
  grandpa.now_init_cond[0] = 0.0;
  grandpa.now_init_cond[1] = 0.02;
  grandpa.now_init_cond[2] = 0.0;
  grandpa.now_init_cond[3] = 0.0;
  grandpa.now_init_cond[4] = 0.0;
  // Perform algorithm
  int i = 1;
  grandpa.check_now_phase();
  while(grandpa.now_init_cond[4] <= 60.0 && i < 60){
    grandpa.init_next_phase(60.0);
    int now_tag = grandpa.now_phase;
    int next_tag = grandpa.next_phase;
    double t_end = grandpa.next_init_cond[4];
    double t_begin = grandpa.now_init_cond[4];
    printf("%1d\t %1d\t %4.7e\t %4.7e\n",now_tag,next_tag,t_begin,t_end);
    i++;
    for(int ii = 0; ii < 5; ii++){
      grandpa.now_init_cond[ii] = grandpa.next_init_cond[ii];
      grandpa.now_phase = grandpa.next_phase;
    }
  }
  // Return memory
  grandpa.kill_punchet();
}
