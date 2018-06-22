#include <cstdlib>
#include <cmath>
#include <fftw3.h>

int main(void){

  int N = 100000;
  double dt = 0.001;
  double w = 2*M_PI;

  fftw_complex *in, *out;
  fftw_plan plan;

  in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N);
  out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N);
  plan = fftw_plan_dft_1d(N,in,out,FFTW_FORWARD,FFTW_ESTIMATE);

  double t = 0.0;
  for(int ii = 0; ii < N; ii++){
    t = ii*dt;
    in[ii][0] = sin(2.0*w*t);
    in[ii][1] = 0.0;
  }

  fftw_execute(plan);

  double f_samp = 1.0/dt;
  double f = 0.0;
  for(int ii = 0; f < 5.0; ii++){
    f = (f_samp*ii)/N;
    double p = out[ii][0]*out[ii][0] + out[ii][1]*out[ii][1];
    printf("%4.5f\t %4.5e\n",f,p);
  }

  fftw_destroy_plan(plan);
  fftw_free(in);
  fftw_free(out);

  return 0;
}
