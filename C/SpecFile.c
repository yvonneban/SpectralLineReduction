#include <stdlib.h>
#include <stdio.h>
#include <netcdf.h>

#include "SpecFile.h"

int read_spec_file(SpecFile *S, char *filename)
{
  int i,j;
  /* This will be the netCDF ID for the file and data variable. */
  int ncid;
  /* for error handling. */
  int retval;
  /* dimensions and varid's */
  size_t nspec, nchan;
  int nspec_id, nchan_id, data_id, x_id, y_id, pix_id, seq_id, rms_id;

  //printf("about to open file %s\n",filename);

  /* Open the file. NC_NOWRITE tells netCDF we want read-only access
     to the file.*/
  if ((retval = nc_open(filename, NC_NOWRITE, &ncid)))
    ERR(retval);

  printf("file %s opened\n",filename);

  /* Get the dimensions */
  if ((retval = nc_inq_dimid(ncid, "nspec", &nspec_id)))
    ERR(retval);
  if ((retval = nc_inq_dimid(ncid, "nchan", &nchan_id)))
    ERR(retval);
  if ((retval = nc_inq_dimlen(ncid, nspec_id, &nspec)))
    ERR(retval);
  if ((retval = nc_inq_dimlen(ncid, nchan_id, &nchan)))
    ERR(retval);

  //printf("Dimensions complete %zu %zu\n",nspec,nchan);

  int obsnum_id, source_id,source_x,source_y,crval_id,crpix_id,cdelt_id,ctype_id,caxis_id;
  /* Get the varids of the observation header */
  if ((retval = nc_inq_varid(ncid, "Header.Obs.ObsNum", &obsnum_id)))
    ERR(retval);
  if ((retval = nc_inq_varid(ncid, "Header.Obs.SourceName", &source_id)))
    ERR(retval);
  if ((retval = nc_inq_varid(ncid, "Header.Obs.XPosition", &source_x)))
    ERR(retval);
  if ((retval = nc_inq_varid(ncid, "Header.Obs.YPosition", &source_y)))
    ERR(retval);

  //printf("Header.Obs complete\n");

  /* Get the varids of the spectrum axis header variables */
  if ((retval = nc_inq_varid(ncid, "Header.SpectrumAxis.CRVAL", &crval_id)))
    ERR(retval);
  //printf("Header.SpectrumAxis.CRVAL\n");
  if ((retval = nc_inq_varid(ncid, "Header.SpectrumAxis.CRPIX", &crpix_id)))
    ERR(retval);
  //printf("Header.SpectrumAxis.CRPIX\n");
  if ((retval = nc_inq_varid(ncid, "Header.SpectrumAxis.CDELT", &cdelt_id)))
    ERR(retval);
  //printf("Header.SpectrumAxis.CDELT\n");
  if ((retval = nc_inq_varid(ncid, "Header.SpectrumAxis.CTYPE", &ctype_id)))
    ERR(retval);
  //printf("Header.SpectrumAxis.CTYPE\n");
  if ((retval = nc_inq_varid(ncid, "Header.SpectrumAxis.CAXIS", &caxis_id)))
    ERR(retval);
  //printf("Header.SpectrumAxis.XAxis\n");

  /* Get the varid of the data variables, based on its name. */
  //printf("Data.Spectra\n");
  if ((retval = nc_inq_varid(ncid, "Data.Spectra", &data_id)))
    ERR(retval);
  //printf("Data.Spectra completed\n");
  if ((retval = nc_inq_varid(ncid, "Data.XPos", &x_id)))
    ERR(retval);
  //printf("Data.XPos\n");
  if ((retval = nc_inq_varid(ncid, "Data.YPos", &y_id)))
    ERR(retval);
  //printf("Data.YPos\n");
  if ((retval = nc_inq_varid(ncid, "Data.Pixel", &pix_id)))
    ERR(retval);
  //printf("Data.Pixel\n");
  if ((retval = nc_inq_varid(ncid, "Data.Sequence", &seq_id)))
    ERR(retval);
  //printf("Data.Sequence\n");
  if ((retval = nc_inq_varid(ncid, "Data.RMS", &rms_id)))
    ERR(retval);
  //printf("Data.RMS\n");

  /* Read the data and load the SpecFile struct */
  printf("file: %s nspec= %zu nchan= %zu\n",filename,nspec,nchan);
  S->nspec = nspec;
  S->nchan = nchan;

  if((retval = nc_get_var_int(ncid,obsnum_id, &S->obsnum)) != NC_NOERR)
    ERR(retval);
  if((retval = nc_get_var(ncid,source_id, S->source)) != NC_NOERR)
    ERR(retval);
  if((retval = nc_get_var_double(ncid, source_x, &S->x_position)) != NC_NOERR)
    ERR(retval);
  if((retval = nc_get_var_double(ncid, source_y, &S->y_position)) != NC_NOERR)
    ERR(retval);

  if((retval = nc_get_var_double(ncid, crval_id, &S->CRVAL)) != NC_NOERR)
    ERR(retval);
  if((retval = nc_get_var_double(ncid, crpix_id, &S->CRPIX)) != NC_NOERR)
    ERR(retval);
  if((retval = nc_get_var_double(ncid, cdelt_id, &S->CDELT)) != NC_NOERR)
    ERR(retval);
  if((retval = nc_get_var(ncid, ctype_id, S->CTYPE)) != NC_NOERR)
    ERR(retval);
  S->CAXIS = (double*)malloc(nchan*sizeof(double));
  if((retval = nc_get_var(ncid, caxis_id, S->CAXIS)) != NC_NOERR)
    ERR(retval);



  S->theData = (float *)malloc(nspec*nchan*sizeof(float));
  if(S->theData == NULL)
    fprintf(stderr,"SpecFile: Error allocating data array\n");
  else
    printf("allocated theData\n");

  S->XPos = (float *)malloc(nspec*sizeof(float));
  S->YPos = (float *)malloc(nspec*sizeof(float));
  S->RMS = (float *)malloc(nspec*sizeof(float));

  S->Pixel = (int *)malloc(nspec*sizeof(int));
  S->Sequence = (int *)malloc(nspec*sizeof(int));

  if ((retval = nc_get_var_float(ncid, data_id, S->theData)))
    ERR(retval);
  if ((retval = nc_get_var_float(ncid, x_id, S->XPos)))
    ERR(retval);
  if ((retval = nc_get_var_float(ncid, y_id, S->YPos)))
    ERR(retval);
  if ((retval = nc_get_var_float(ncid, rms_id, S->RMS)))
    ERR(retval);
  if ((retval = nc_get_var_int(ncid, pix_id, S->Pixel)))
    ERR(retval);
  if ((retval = nc_get_var_int(ncid, seq_id, S->Sequence)))
    ERR(retval);

  /* Close the file, freeing all resources. */
  if ((retval = nc_close(ncid)))
    ERR(retval);
  return 0;
}

void free_spec_file(SpecFile *S)
{
  free(S->Sequence);
  free(S->Pixel);
  free(S->RMS);
  free(S->YPos);
  free(S->XPos);
  free(S->theData);
  free(S->CAXIS);
}

float *get_spectrum(SpecFile *S, int i)
{
  int index;

  index = i*S->nchan;
  return(&(S->theData[index]));
}
