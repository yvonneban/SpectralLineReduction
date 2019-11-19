#ifndef _CONVOLVE_H_
#define _CONVOLVE_H_

typedef struct
{
  int nc,ns,nn,xc,yc,nsub;
  float cell_size;
  float subcell_size;
  float **array;
} ConvolveArray;

void initialize_convolve_array(ConvolveArray *CA, int nc, int ns, float cell_size);
void initialize_box(ConvolveArray *CA);
void initialize_gaussian(ConvolveArray *CA, float hpbw);
void initialize_jinc(ConvolveArray *CA, float, float, float);
void print_subcell_grid(ConvolveArray *CA);
void print_convolve_sub_array(ConvolveArray *CA, int istart, int jstart);
void print_convolve_array(ConvolveArray *CA);
#endif
