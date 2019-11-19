#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "ConvolveArray.h"
double jinc(double r);

void initialize_convolve_array(ConvolveArray *CA, int nc, int ns, float cell_size)
{
  int i,j;
  CA->nc = nc;
  CA->ns = ns;
  CA->nn = (2*CA->nc+1)*(2*CA->ns+1);
  CA->nsub = 2*CA->ns+1;
  CA->cell_size = cell_size;
  CA->subcell_size = CA->cell_size / CA->nsub;

  CA->xc = (CA->nn-1)/2;
  CA->yc = (CA->nn-1)/2;

  // initialize convolution array with zeros to start
  CA->array = (float**)malloc(CA->nn*sizeof(float*));
  for(i=0;i<CA->nn;i++)
    {
      CA->array[i] = (float*)malloc(CA->nn*sizeof(float));
      for(j=0;j<CA->nn;j++)
	CA->array[i][j] = 0.0;
    }
}

void print_subcell_grid(ConvolveArray *CA)
{
  int i,j;
  float x,y;
  for(i=0; i<CA->nn; i++)
    {
      x = i*CA->subcell_size - CA->xc*CA->subcell_size;
      for(j=0; j<CA->nn; j++)
	{
	  y = j*CA->subcell_size - CA->yc*CA->subcell_size;
	  printf("%d %d %f %f\n",i,j,x,y);
	}
    }
}

void print_convolve_sub_array(ConvolveArray *CA, int istart, int jstart)
{
  int i,j;

  printf("   ");
  for(j=jstart;j<CA->nn; j+=CA->nsub)
    printf("%5d ",j);
  printf("\n");
  for(i=istart;i<CA->nn;i+=CA->nsub)
    {
      printf("%2d ",i);
      for(j=jstart;j<CA->nn;j+=CA->nsub)
	{
	  printf("%5.1f ",CA->array[i][j]);
	}
      printf("\n");
    }
}

void print_convolve_array(ConvolveArray *CA)
{
  int i,j,ii,jj;
  int idx,idy;
  float xx,yy;
  for(i=-CA->nc; i<=CA->nc; i++)
    for(j=-CA->nc; j<=CA->nc; j++)
      for(ii=-CA->ns; ii<=CA->ns; ii++)
	for(jj=-CA->ns; jj<=CA->ns; jj++)
	  {
	    xx = i*CA->cell_size - ii*CA->subcell_size;
	    yy = j*CA->cell_size - jj*CA->subcell_size;
	    idx = (i+CA->nc)*(2*CA->ns+1) + ii + CA->ns;
	    idy = (j+CA->nc)*(2*CA->ns+1) + jj + CA->ns;
	    printf("%2d %2d %2d %2d  %2d %2d   %5.1f %5.1f   %6.3f\n",i,j,ii,jj,idx,idy,xx,yy,CA->array[idx][idy]);
	    
	  }

}

void initialize_box(ConvolveArray *CA)
{
  int i,j;

  // center pixel
  for(i=CA->xc-CA->ns; i<=CA->xc+CA->ns; i++)
    for(j=CA->yc-CA->ns; j<=CA->yc+CA->ns; j++)
      CA->array[i][j] = 1.0;
}

void initialize_gaussian(ConvolveArray *CA, float hp_width)
{
  int i,j;
  int n;
  float x,y, xx,yy;
  int ii,jj,idx,idy;

  for(i=-CA->nc; i<=CA->nc; i++)
    for(j=-CA->nc; j<=CA->nc; j++)
      for(ii=-CA->ns; ii<=CA->ns; ii++)
	for(jj=-CA->ns; jj<=CA->ns; jj++)
	  {
	    xx = i*CA->cell_size - ii*CA->subcell_size;
	    yy = j*CA->cell_size - jj*CA->subcell_size;
	    idx = (i+CA->nc)*(2*CA->ns+1) + ii + CA->ns;
	    idy = (j+CA->nc)*(2*CA->ns+1) + jj + CA->ns;
	    CA->array[idx][idy] = exp(-2.7726*(xx*xx+yy*yy)/(hp_width*hp_width));	    
	  }
}

// we use the function from Mangum, Emerson, and Greisen A&A 474 679-687 (2007)
// nominally the constants are:
// otf_a = 1.55 * hpbw/3.
// otf_b = 2.52 * hpbw/3.
// otf_c = 2.

void initialize_jinc(ConvolveArray *CA, float a_otf, float b_otf, float c_otf)
{
  int i,j;
  int n;
  float x,y, xx, yy, rr;
  int ii,jj,idx,idy;

  for(i=-CA->nc; i<=CA->nc; i++)
    for(j=-CA->nc; j<=CA->nc; j++)
      for(ii=-CA->ns; ii<=CA->ns; ii++)
	for(jj=-CA->ns; jj<=CA->ns; jj++)
	  {
	    xx = i*CA->cell_size - ii*CA->subcell_size;
	    yy = j*CA->cell_size - jj*CA->subcell_size;
	    idx = (i+CA->nc)*(2*CA->ns+1) + ii + CA->ns;
	    idy = (j+CA->nc)*(2*CA->ns+1) + jj + CA->ns;
	    rr = sqrt(xx*xx+yy*yy);
	    if(rr == 0.0)
	      CA->array[idx][idy] = exp(-pow(rr/b_otf,c_otf));
	    else
	      CA->array[idx][idy] = jinc(rr/a_otf)/(rr/a_otf)*exp(-pow((rr/b_otf),c_otf));
	  }
}



// jinc function taken from the fcrao otf program
double jinc(double r)
{
  double arg,f1,theta,j;

  if (r <=3.0) 
    {
      arg=pow(r/3.0,2.0);
      j=0.5+arg*(-0.56249985+arg*(0.21093573+arg*(-0.03954289+arg*(0.00443319+arg*(-0.00031761+arg*0.00001109)))));
    }
  else 
    {
      arg=3.0/r;
      f1=0.79788456+arg*(0.00000156+arg*(0.01659667+arg*(0.00017105+arg*(-.00249511+arg*(0.00113653-arg*0.00020033)))));
      theta=r-2.35619449+arg*(0.12499612+arg*(0.00005650+arg*(-.00637879+arg*(0.00074348+arg*(0.00079824-arg*0.00029166))))) ;
   
      j=f1*cos(theta)/(r*sqrt(r));
    }
  return(j);
}
