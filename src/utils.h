#ifndef UTILS_H
#define UTILS_H

#include <stdlib.h>


size_t read_video(
  char *name,       //file name
  unsigned char *I, //video to read
  int size          //size of the video
);

size_t write_video(
  char *name,       //file name
  unsigned char *I, //video to write
  int size          //size of the video
);


void save_transforms(
  char *name,      //file name
  float *H,        //set of transformations
  int nparams,     //number of parameters of the transformations
  int ntransforms, //number of transformations
  int nx,          //image width
  int ny           //image height
);

void read_transforms(
  char *name,      //file name
  float *H,        //set of transformations
  int nparams,     //number of parameters of the transformations
  int ntransforms, //number of transformations
  int &nx,         //image width
  int &ny          //image height
);


#endif
