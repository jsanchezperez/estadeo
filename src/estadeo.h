#ifndef ESTADEO_H
#define ESTADEO_H


/**
  *
  * Function for Online Video Estabilization
  * 
**/
void estadeo_online(
  float *I,            //input grayscale video to estabilize
  float *Ic,           //input color video to estabilize
  int nx,              //number of columns 
  int ny,              //number of rows
  int nz,              //number of channels
  int nframes,         //number of frames of the video
  int nparams,         //type of matrix transformation
  int smooth_strategy, //motion smoothing strategy
  float *sigma,        //Gaussian standard deviations
  char *out_stransform,//output file to write the stabilizing transformations
  int verbose          //verbose mode
);

#endif
