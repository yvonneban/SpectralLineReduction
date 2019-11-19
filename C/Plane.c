/** plane.c - methods for handling data plane
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "Plane.h"

/*************************  PLANE METHODS ***********************************/

/** initialize_plane - allocated data for the plane
 */
void initialize_plane(Plane* P, int *n)
{
  int i;
  for(i=0;i<2;i++)
    P->n[i] = n[i];
  P->nplane = n[0]*n[1];
  P->plane = (float*)malloc(P->nplane*sizeof(float));
  if(P->plane == NULL)
    fprintf(stderr,"Plane: Failed to Allocate Plane\n");
  else
    for(i=0;i<P->nplane;i++)
      P->plane[i] = 0.0;
}

/** initialize_axis - initializes the axis data 
 */
void initialize_plane_axis(Plane *P, int axis, float crval, float crpix, float cdelt, char *ctype, char *cunit)
{
  int i;

  P->crval[axis] = crval;
  P->crpix[axis] = crpix;
  P->cdelt[axis] = cdelt;
  strcpy(P->ctype[axis],ctype);
  strcpy(P->cunit[axis],cunit);

  P->caxis[axis] = (float*)malloc(P->n[axis]*sizeof(float));
  if(P->caxis[axis] == NULL)
    fprintf(stderr,"Plane: Failed to allocate axis %d\n",axis);

  for(i=0;i<P->n[axis];i++)
      P->caxis[axis][i] = (i-P->crpix[axis])*P->cdelt[axis]+P->crval[axis];
}

/** axis_index - finds the index in an axis array given a value
 */
int plane_axis_index(Plane *P, int axis, float value)
{
  float result_f;
  int result_i;
  
  result_f = (value-P->crval[axis])/P->cdelt[axis] + P->crpix[axis];
  result_i = (int)floor(result_f);
  if((result_i<0) || (result_i>=P->n[axis]))
    result_i = -1;
  return(result_i);
}

/** plane_index - finds the index of an element in the plane according to 
    given x and y positions
*/
int plane_index(Plane *P, float x, float y)
{
  int ix, iy;
  int result;

  ix = plane_axis_index(P, PLANE_X_AXIS, x);
  iy = plane_axis_index(P, PLANE_Y_AXIS, y);

  result = iy*P->n[PLANE_X_AXIS] + ix;
  return(result);
}

