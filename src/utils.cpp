// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2017, Javier Sánchez Pérez <jsanchez@ulpgc.es>
// All rights reserved.


#include <stdio.h>


/**
  *
  *  Function to read a video in raw data
  * 
**/
size_t read_video(
  char *name,       //file name
  unsigned char *I, //video to read
  int size          //size of the video
) 
{
  FILE* file=fopen(name, "rb");
  if(file != NULL)
  {
    size_t result=fread(I, sizeof(unsigned char), size, file);
    fclose(file);
    return result;
  }
  else return 0;
}


/**
  *
  *  Function to read a video in raw data
  * 
**/
size_t write_video(
  char *name,       //file name
  unsigned char *I, //video to write
  int size          //size of the video
) 
{
  FILE* file=fopen(name, "wb");
  if(file != NULL)
  {
    size_t result=fwrite(I, sizeof(unsigned char), size, file);
    fclose(file);
    return result;
  }
  else return 0;
}


/**
  *
  *  Function to save transformations to a file
  * 
**/
void save_transforms(
  char *name,      //file name
  float *H,        //set of transformations
  int nparams,     //number of parameters of the transformations
  int ntransforms, //number of transformations
  int nx,          //image width
  int ny           //image height
)
{
  FILE *fd=fopen(name,"w");
  if(fd!=NULL)
  {
    fprintf(fd,"%d %d %d %d\n", nparams, ntransforms, nx, ny);

    for(int i=0;i<ntransforms;i++)
    {
      for(int j=0;j<nparams;j++) fprintf(fd, "%.15lf ", H[i*nparams+j]);
      fprintf(fd, "\n");
    }
    fclose(fd);
  }
}


/**
  *
  *  Function to read transformations from a file
  * 
**/
void read_transforms(
  char *name,      //file name
  float *H,        //set of transformations
  int nparams,     //number of parameters of the transformations
  int ntransforms, //number of transformations
  int &nx,         //image width
  int &ny          //image height
)
{
  FILE *fd=fopen(name,"r");
  if(fd!=NULL)
  {
    int r=fscanf(fd,"%d %d %d %d", &nparams, &ntransforms, &nx, &ny);
    if(r>0)
      for(int i=0;i<ntransforms;i++) {
	for(int j=0;j<nparams;j++) r=fscanf(fd, "%f", &(H[i*nparams+j]));
      }
    fclose(fd);
  }
}
