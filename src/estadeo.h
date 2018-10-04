#ifndef ESTADEO_H
#define ESTADEO_H

#define DIRECT_METHOD 0
#define FEATURE_BASED_METHOD 1


/**
  *
  * Function for Offline Video Estabilization
  *
**/
void estadeo(
  float *I,            //input grayscale video to estabilize
  float *Ic,           //input color video to estabilize
  int nx,              //number of columns 
  int ny,              //number of rows
  int nz,              //number of channels
  int nframes,         //number of frames of the video
  int nparams,         //type of matrix transformation
  int smooth_strategy, //motion smoothing strategy
  float *sigma,        //Gaussian standard deviations
  int bc,              //boundary condition for smoothing
  int postprocessing,  //method for dealing with empty regions
  char *in_transform,  //input file with the computed transformations
  char *out_transform, //output file to write the computed transformations
  char *out_stransform,//output file to write the stabilizing transformations
  int verbose          //verbose mode
);


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
  int bc,              //boundary condition for smoothing
  int postprocessing,  //method for dealing with empty regions
  char *out_stransform,//output file to write the stabilizing transformations
  int verbose          //verbose mode
);

#endif
