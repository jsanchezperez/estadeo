// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2017-2018, Javier Sánchez Pérez <jsanchez@ulpgc.es>
// All rights reserved.


#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h> 

#include "estadeo.h"
#include "motion_smoothing.h"
#include "utils.h"
#include "color_bicubic_interpolation.h"
#include "inverse_compositional_algorithm.h"
#include "transformation.h"



/**
  *
  * Function for estimating the transformation between two frames
  *
**/
void online_motion
(
  float *I1,   //first image
  float *I2,   //second image
  float *H,    //output matrix transformation
  int nparams, //type of matrix transformation
  int   nx,    //number of columns
  int   ny     //number of rows
)
{
  //parameters for the direct method
  int   nscales=10;
  float TOL=1E-3;
  float lambda=0;
  float N;
  int   robust=LORENTZIAN; //QUADRATIC; //
  int   max_d=200;
  int   min_d=50;

  N=1+log(((nx<ny)?nx:ny)/min_d)/log(2.);
  if ((int) N<nscales) nscales=(int) N;

  //motion estimation through direct methods
  pyramidal_inverse_compositional_algorithm(
    I1, I2, H, nparams, nx, ny, nscales, TOL, robust, lambda, max_d
  );
}


/**
  *
  * Function for onlie warping the last frame of the video
  *
**/
void online_warp
(
  float *I,    //frame to be warped
  float *H,    //matrix transformations
  int nparams, //type of matrix transformation
  int nx,      //number of columns   
  int ny,      //number of rows
  int nz       //number of channels
)
{
  int size=nx*ny*nz;

  float *I2=new float[nx*ny*nz];

  //warp the image
  //bicubic_interpolation
  bilinear_interpolation(I, I2, H, nparams, nx, ny, nz);

  //copy warped image
  for(int j=0; j<size; j++)
    I[j]=I2[j];

  delete []I2;
}


/**
  *
  * Function for Online Video Estabilization
  * 
**/
void estadeo_online(
  float *I,       //input grayscale video to estabilize
  float *Ic,      //input color video to estabilize
  int   nx,       //number of columns 
  int   ny,       //number of rows
  int   nz,       //number of channels
  int   nframes,  //number of frames of the video
  int   nparams,  //type of matrix transformation
  int   strategy, //motion smoothing strategy
  float sigma,    //Gaussian standard deviation
  char  *Hout,    //output file to write the stabilizing transformations
  int   verbose   //verbose mode
)
{ 
  int size=nx*ny;
  float *H=new float[nframes*nparams];
  float *Hp=new float[nframes*nparams];
  float avg1=0, avg2=0, avg3=0;

  //introduce identity matrix for the first transform
  for(int i=0; i<nparams; i++) H[i]=0;

  for(int f=1; f<nframes; f++)
  {
    struct timeval t1, t2, t3, t4;
    
    //step 1. Compute motion between the last two frames
    gettimeofday(&t1, NULL);
    online_motion(
      &I[size*(f-1)], &I[size*f], &H[nparams*f], nparams, nx, ny
    );
    
    //step 2. Smooth until the last transformation 
    gettimeofday(&t2, NULL);
    online_smoothing(H, Hp, nparams, f+1, sigma, strategy);
    
    //step 3. Warp the image
    gettimeofday(&t3, NULL);
    online_warp(&Ic[size*nz*f], &Hp[nparams*f], nparams, nx, ny, nz);
    
    if(verbose){
      gettimeofday(&t4, NULL);
      avg1+=((t2.tv_sec-t1.tv_sec)*1000000u+t2.tv_usec-t1.tv_usec)/1.e6;
      avg2+=((t3.tv_sec-t2.tv_sec)*1000000u+t3.tv_usec-t2.tv_usec)/1.e6;
      avg3+=((t4.tv_sec-t3.tv_sec)*1000000u+t4.tv_usec-t3.tv_usec)/1.e6;
      printf(" Processing frame %d: T(%.4fs, %.7fs, %.4fs) \n", f, 
            ((t2.tv_sec-t1.tv_sec)*1000000u+t2.tv_usec-t1.tv_usec)/1.e6,
            ((t3.tv_sec-t2.tv_sec)*1000000u+t3.tv_usec-t2.tv_usec)/1.e6,
            ((t4.tv_sec-t3.tv_sec)*1000000u+t4.tv_usec-t3.tv_usec)/1.e6);
    }
  }

  if(verbose){
    float total=avg1+avg2+avg3;
    float average=(avg1+avg2+avg3)/nframes;
    printf(
      "\n Average time per frame: %.4fs -> T(%.4fs, %.7fs, %.4fs) \n", 
      average, avg1/nframes, avg2/nframes, avg3/nframes
    ); 
    printf(
      "\n Total time: %.4fs -> T(%.4fs, %.7fs, %.4fs) \n", 
      total, avg1, avg2, avg3
    ); 
  }
  
  //save the stabilizing transformations 
  if(out_stransform!=NULL)
  {
    if(verbose)
      printf(
        "\n  Write smoothed transformations to file '%s'\n", out_stransform
      );

    float *Ho=new float[nframes*nparams];
    compute_smooth_transforms(H, Hp, Ho, nparams, nframes);
    save_transforms(out_stransform, Ho, nparams, nframes, nx, ny);
    delete []Ho;
  }

  delete []H;
  delete []Hp;
}


