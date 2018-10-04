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
#include "crop_and_zoom.h"
#include "motion_smoothing.h"
#include "utils.h"
#include "color_bicubic_interpolation.h"
#include "inverse_compositional_algorithm.h"
#include "transformation.h"


#define ICA //ICA_2

/**
  *
  * Function for estimating the transformation between frames
  *
**/
void motion_estimation
(
  float *I,        //input video to estabilize
  float *H,        //output matrix transformations
  int nparams,     //type of matrix transformation
  int nx,          //number of columns
  int ny,          //number of rows
  int nz,          //number of channels
  int nframes,     //number of frames of the video
  int verbose      //verbose mode
)
{
  //introduce identity matrix for the first transform
  for(int i=0;i<nparams;i++) H[i]=0;

  int size=nx*ny*nz;

  if(verbose) printf("  Direct method: ");
 
  //parameters for the direct method
  int nscales=10;
  float zfactor=0.5;
  float TOL=1E-3;
  float lambda=0;
  float N;
  int robust=QUADRATIC;  //LORENTZIAN;

  N=1+log(((nx<ny)?nx:ny)/50.)/log(1./zfactor);
  if ((int) N<nscales) nscales=(int) N;

  //motion estimation through direct methods
  #pragma omp parallel for
  for(int i=0;i<nframes-1;i++) 
  {
    if(verbose) 
    {
      printf("%d, ",i);
      fflush(stdout);
    } 

    pyramidal_inverse_compositional_algorithm(
      &I[size*i], &I[size*(i+1)], &(H[(i+1)*nparams]), nparams, 
      nx, ny, nz, nscales, zfactor, TOL, robust, lambda
    );
  }

  if(verbose) printf("\n");
}


/**
  *
  * Function for warping the frames of the input video
  *
**/
void warp_video
(
  float *I,    //video to be warped
  float *H,    //matrix transformations
  int nparams, //type of matrix transformation
  int nx,      //number of columns   
  int ny,      //number of rows
  int nz,      //number of channels
  int nframes, //number of frames of the video
  int border,  //do not fill the border
  int verbose  //verbose mode
)
{
  int size=nx*ny*nz;

  if(verbose)
    printf("  ");

  //warp the frames
  #pragma omp parallel for
  for(int i=0;i<nframes;i++)
  {
    if(verbose)
    {
      printf("%d, ", i);
      fflush(stdout);
    }

    float *I2=new float[nx*ny*nz];

    //warp the image
    bicubic_interpolation(
      &(I[size*i]), I2, &(H[i*nparams]), nparams, nx, ny, nz, border==0
    );

    //copy warped image
    for(int j=0; j<size; j++)
      I[size*i+j]=I2[j];

    delete []I2;
  }

  if(verbose) printf("\n");
}


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
  int nparams,  //type of matrix transformation
  int nx,       //number of columns
  int ny,        //number of rows
  int nz        //number of channels
)
{
  //parameters for the direct method
  int nscales=10;
  float zfactor=0.5;
  float TOL=1E-3;
  float lambda=0;
  float N;
  int robust=QUADRATIC; //LORENTZIAN;

  N=1+log(((nx<ny)?nx:ny)/50.)/log(1./zfactor);
  if ((int) N<nscales) nscales=(int) N;

  //motion estimation through direct methods
  pyramidal_inverse_compositional_algorithm(
    I1, I2, H, nparams, nx, ny, nz, nscales, zfactor, TOL, robust, lambda
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
  int nz,      //number of channels
  int border   //do not fill the border
)
{
  int size=nx*ny*nz;

  float *I2=new float[nx*ny*nz];

  //warp the image
  bicubic_interpolation(I, I2, H, nparams, nx, ny, nz, border==0);

  //copy warped image
  for(int j=0; j<size; j++)
    I[j]=I2[j];

  delete []I2;
}


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
)
{
  float *H=new float[nframes*nparams];
  float *Hp=new float[nframes*nparams];

  //--------------------------------
  //motion estimation between frames
  //--------------------------------
  if(in_transform==NULL)
  {
    if(verbose)
      printf("\n 1-Motion estimation:\n");

    //compute the motion between consecutive frames
#ifdef ICA  //use ICA for graylevel images
    motion_estimation(I, H, nparams, nx, ny, 1, nframes, verbose);
#else
    motion_estimation(Ic, H, nparams, nx, ny, nz, nframes, verbose);
#endif

    //save the original transformation
    if(out_transform!=NULL)
    {
      if(verbose)
        printf("\n  Write transformations to file '%s'\n", out_transform);

      save_transforms(out_transform, H, nparams, nframes, nx, ny);
    }
  }
  else
  {
    //load transformations from the input file
    if(verbose)
      printf(" 1-Loading motion parameters from '%s':\n", in_transform);

    read_transforms(in_transform, H, nparams, nframes, nx ,ny);
  }

  //----------------------------------------------------------------  
  //smooth the transformations and obtain the stabilizing parameters
  //----------------------------------------------------------------
  if(verbose)
    printf("\n 2-Motion smoothing:\n");

  motion_smoothing(
    H, Hp, nparams, nframes, sigma, smooth_strategy, bc, verbose
  );

  //save the stabilizing transformation 
  if(out_stransform!=NULL)
  {
    if(verbose)
      printf("\n  Write smoothed transformations to file '%s'\n", 
             out_stransform);

    float *Ho=new float[nframes*nparams];
    compute_smooth_transforms(H, Hp, Ho, nparams, nframes);
    save_transforms(out_stransform, Ho, nparams, nframes, nx, ny);
    delete []Ho;
  }


  //---------------------------------------------
  //postprocessing step for removing empty spaces
  //---------------------------------------------
  if(verbose)
    printf("\n 3-Postprocessing:\n");

  crop_and_zoom(Hp, nx, ny, nframes, nparams, postprocessing, verbose);


  //-------------------------------------------------------
  //warp the input video to produce the stabilized sequence
  //-------------------------------------------------------
  if(verbose)
    printf("\n 4-Video warping:\n");

  warp_video(
    Ic, Hp, nparams, nx, ny, nz, 
    nframes, postprocessing, verbose
  );

  delete []H;
  delete []Hp;
}



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
)
{
  int size=nx*ny;
  float *H=new float[nframes*nparams];
  float *Hp=new float[nframes*nparams];

  //introduce identity matrix for the first transform
  for(int i=0;i<nparams;i++) H[i]=0;

  for(int f=1; f<nframes; f++)
  {
    struct timeval t1, t2, t3, t4;
    
    //step 1. Compute motion between the last two frames
    gettimeofday(&t1, NULL);
#ifdef ICA  //use ICA for graylevel images
    online_motion(
      &I[size*(f-1)], &I[size*f], &H[nparams*f], nparams, nx, ny, 1
    );
#else  //no está funcionando bien.........
    online_motion(
      &Ic[size*nz*(f-1)], &Ic[size*nz*f], &H[nparams*f], nparams, nx, ny, nz
    );
#endif 
    
    //step 2. Smooth the last transformation 
    gettimeofday(&t2, NULL);
    online_smoothing(H, Hp, nparams, f+1, sigma, smooth_strategy, bc);
    
    //step 3. Warp the image
    gettimeofday(&t3, NULL);
    online_warp(
      &Ic[size*nz*f], &Hp[nparams*f], nparams, nx, ny, nz, postprocessing
    );
    
    if(verbose){
      gettimeofday(&t4, NULL);
      printf(" Processing frame %d: T(%.4fs, %.4fs, %.4fs) \n", f, 
            ((t2.tv_sec-t1.tv_sec)*1000000u+t2.tv_usec-t1.tv_usec)/1.e6,
            ((t3.tv_sec-t2.tv_sec)*1000000u+t3.tv_usec-t2.tv_usec)/1.e6,
            ((t4.tv_sec-t3.tv_sec)*1000000u+t4.tv_usec-t3.tv_usec)/1.e6);
    }
  }

  //save the stabilizing transformations 
  if(out_stransform!=NULL)
  {
    if(verbose)
      printf("\n  Write smoothed transformations to file '%s'\n", 
             out_stransform);

    float *Ho=new float[nframes*nparams];
    compute_smooth_transforms(H, Hp, Ho, nparams, nframes);
    save_transforms(out_stransform, Ho, nparams, nframes, nx, ny);
    delete []Ho;
  }

  delete []H;
  delete []Hp;
}


