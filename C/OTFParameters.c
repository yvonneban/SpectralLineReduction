#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <assert.h>
#include "OTFParameters.h"

char** str_split(char* a_str, const char a_delim)
{
    char** result    = 0;
    size_t count     = 0;
    char* tmp        = a_str;
    char* last_comma = 0;
    char delim[2];
    delim[0] = a_delim;
    delim[1] = 0;

    /* Count how many elements will be extracted. */
    while (*tmp)
    {
        if (a_delim == *tmp)
        {
            count++;
            last_comma = tmp;
        }
        tmp++;
    }

    /* Add space for trailing token. */
    count += last_comma < (a_str + strlen(a_str) - 1);

    /* Add space for terminating null string so caller
       knows where the list of returned strings ends. */
    count++;

    result = malloc(sizeof(char*) * count);

    if (result)
    {
        size_t idx  = 0;
        char* token = strtok(a_str, delim);

        while (token)
        {
            assert(idx < count);
            *(result + idx++) = strdup(token);
            token = strtok(0, delim);
        }
        assert(idx == count - 1);
        *(result + idx) = 0;
    }

    return result;
}

void decode_pix_list(OTFParameters *OTF, char *the_list)
{
    int pixel,i,len;
    char p[40];
    char** tokens;

    //printf("pix_list=%s\n\n", the_list);

    len = strlen(the_list);
    strncpy(p,&the_list[1],len-2);
    p[len-1] = '\0';
    tokens = str_split(p, ',');

    if (tokens)
    {
        for (i = 0; *(tokens + i); i++)
        {
            pixel = atoi(*(tokens + i));
            OTF->use_pixels[pixel] = 1;
            //printf("pixel=%d\n", pixel);
            free(*(tokens + i));
        }
	// printf("\n");
        free(tokens);
    }
    /*
    for(i=0;i<16;i++)
	printf("pixel %d use_it %d\n",i,OTF->use_pixels[i]);
    */
}

void decode_file_list(OTFParameters *OTF, char *the_list)
{
    int i;
    char** tokens;

    printf("file list=%s\n", the_list);

    tokens = str_split(the_list, ',');

    OTF->nfiles = 0;
    if (tokens)
    {
        for (i = 0; *(tokens + i); i++)
        {
	  OTF->nfiles++;
          strcpy(OTF->i_filename[i],*(tokens + i));
	  free(*(tokens + i));
        }
        free(tokens);
    }
}


void initialize_otf_parameters(OTFParameters *OTF, int argc, char *argv[])
{
  int i;
  int coption;
  /*
    Default Parameters
   */
  strcpy(OTF->o_filename,"test.nc");
  OTF->rms_cutoff = 10.;
  OTF->noise_sigma = 0.0;
  OTF->resolution_size = 14.;
  OTF->cell_size = 7.;
  OTF->nsamples = 256;
  OTF->n_cell = 5;
  OTF->n_subcell = 2;
  OTF->model_spectrum_hpw = 5.;
  OTF->model_source_amp = 1.;
  OTF->model_source_hpw = 16.1;
  OTF->model_source_x = 0.;
  OTF->model_source_y = 0.;
  OTF->model_nchan = 64;

  OTF->otf_jinc_a = 1.1;
  OTF->otf_jinc_b = 4.75;
  OTF->otf_jinc_c = 2.0;
  OTF->otf_select = 1; // 1=jinc other=gauss
  OTF->rmax = 3.0;

  OTF->x_extent = 160.0;
  OTF->y_extent = 160.0;
  OTF->sample_step = 1.0;
  OTF->scan_step = 6.65;
  for(i=0;i<16;i++)
    OTF->use_pixels[i] = 0;
 
  // parse command line arguments
  opterr = 0;

  while(1)
    {
      static struct option long_options[] =
	{
	  {"help",  no_argument, 0, 'h'},
	  {"input",  required_argument, 0, 'i'},
	  {"output",  required_argument, 0, 'o'},
	  {"resolution_size",  required_argument, 0, 'l'},
	  {"cell_size",  required_argument, 0, 'c'},
	  {"rms_cutoff", required_argument, 0, 'z'},
	  {"filter",  required_argument, 0, 'f'},
	  {"n_cell",  required_argument, 0, 'n'},
	  {"n_subcell",  required_argument, 0, 'm'},
	  {"rmax", required_argument, 0, 'r'},
	  {"filter", required_argument, 0, 'f'},
	  {"jinc_a",  required_argument, 0, '0'},
	  {"jinc_b",  required_argument, 0, '1'},
	  {"jinc_c",  required_argument, 0, '2'},
	  {"noise_sigma",  required_argument, 0, 's'},
	  {"x_extent",  required_argument, 0, 'x'},
	  {"y_extent",  required_argument, 0, 'y'},
	  {"sample_step",  required_argument, 0, 'p'},
	  {"scan_step",  required_argument, 0, 'q'},
	  {"pix_list", required_argument, 0, 'u'},
	  {0,0,0,0}
	};

      int option_index=0;
      coption = getopt_long(argc, argv, "hi:o:bjl:z:c:n:f:m:s:r:0:1:2:x:y:p:q:u:",long_options,&option_index);
      
      if(coption == -1)
	break;

      printf("%c\n",coption);
      switch(coption)
	{
	case 'h':
	  for(i=0;i<14;i++)
	    printf("%c %s\n",long_options[i].val,long_options[i].name);
	  exit(0);
	case 'i':
	  decode_file_list(OTF, optarg);
	  //strcpy(OTF->i_filename,optarg);
	  break;
	case 'o':
	  strcpy(OTF->o_filename,optarg);
	  break;
	case 'f':
	  OTF->otf_select = atoi(optarg);
	  break;
	case 'l':
	  OTF->resolution_size = atof(optarg);
	  break;
	case 'c':
	  OTF->cell_size = atof(optarg);
	  break;
	case 'u':
	  decode_pix_list(OTF,optarg);
	  break;
	case 'z':
	  OTF->rms_cutoff = atof(optarg);
	  break;
	case 'n':
	  OTF->n_cell = atoi(optarg);
	  break;
	case 'm':
	  OTF->n_subcell = atoi(optarg);
	  break;
	case 's':
	  OTF->noise_sigma = atof(optarg);
	  break;
	case 'r':
	  OTF->rmax = atof(optarg);
	  break;
	case '0':
	  OTF->otf_jinc_a = atof(optarg);
	  break;
	case '1':
	  OTF->otf_jinc_b = atof(optarg);
	  break;
	case '2':
	  OTF->otf_jinc_c = atof(optarg);
	  break;
	case 'x':
	  OTF->x_extent = atof(optarg);
	  break;
	case 'y':
	  OTF->y_extent = atof(optarg);
	  break;
	case 'p':
	  OTF->sample_step = atof(optarg);
	  break;
	case 'q':
	  OTF->scan_step = atof(optarg);
	  break;
	default:
	  printf("Command Line Argument Problem\n");
	  abort();
	}
    }
  // Derived Values
  OTF->nx_samples = (int)(OTF->x_extent/OTF->sample_step);
  OTF->ny_samples = (int)(OTF->y_extent/OTF->scan_step);
}
