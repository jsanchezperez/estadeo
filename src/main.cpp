// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2017-2018, Javier Sánchez Pérez <jsanchez@ulpgc.es>
// All rights reserved.

#include "estadeo.h"
#include "motion_smoothing.h"
#include "crop_and_zoom.h"
#include "utils.h"
#include "transformation.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define PAR_DEFAULT_OUTVIDEO "output_video.raw"
#define PAR_DEFAULT_TRANSFORM SIMILARITY_TRANSFORM
#define PAR_DEFAULT_SMOOTHING LOCAL_MATRIX_BASED_SMOOTHING
#define PAR_DEFAULT_BILATERAL BILATERAL_COMPOSED
#define PAR_DEFAULT_SIGMA_T 30.0
#define PAR_DEFAULT_SIGMA_H 100.0
#define PAR_DEFAULT_BC NEUMANN_BC
#define PAR_DEFAULT_POSTPROCESS NO_POSTPROCESS
#define PAR_DEFAULT_OUTTRANSFORM "transform.mat"
#define PAR_DEFAULT_INTRANSFORM "in_transform.mat"
#define PAR_DEFAULT_VERBOSE 0


/**
 *
 *  Print a help message 
 *
 */
void print_help(char *name)
{
  printf("\n  Usage: %s raw_input_video width height nframes [OPTIONS] \n\n",
          name);
  printf("  Video stabilization:\n");
  printf("  'raw_input_video' is a video file in raw format (rgb24).\n");
  printf("  'width' is the width of the images in pixels.\n");
  printf("  'height' is the height of the images in pixels.\n");
  printf("  'nframes' is the number of frames in the video.\n");
  printf("  -----------------------------------------------\n");
  printf("  Converting to raw data:\n");
  printf("  'avconv -i video.mp4 -f rawvideo -pix_fmt rgb24 -y "
         "raw_video.raw'\n");
  printf("  to convert an mp4 video to raw format.\n");
  printf("  'avconv -f rawvideo -pix_fmt rgb24 -video_size 640x360 "
         "-framerate\n");
  printf("  30 -i output_video.raw -pix_fmt yuv420p output_video.mp4'\n");
  printf("  to convert a raw video to mp4 format.\n");
  printf("  -----------------------------------------------\n");
  printf("  More information in http://www.ipol.im \n\n");
  printf("  OPTIONS:\n"); 
  printf("  --------\n");
  printf("   -o name  output video name to write the computed raw video\n");
  printf("              default value '%s'\n", PAR_DEFAULT_OUTVIDEO);
  printf("   -t N     transformation type to be computed:\n");
  printf("              2.translation; 3.Euclidean transform;\n");
  printf("              4.similarity; 6.affinity; 8.homography\n"); 
  printf("              default value %d\n", PAR_DEFAULT_TRANSFORM);
  printf("   -m N     motion smoothing strategy:\n");
  printf("              0.pure composition;\n");
  printf("              1.compositional smoothing;\n");
  printf("              2.compositional local smoothing; \n");
  printf("              3.local matrix-based smoothing;\n");
  printf("              4.local point-based smoothing\n");
  printf("              5.local linear matrix-based smoothing\n");
  printf("              6.local linear point-based smoothing\n");
  printf("              default value %d\n", PAR_DEFAULT_SMOOTHING);
  printf("   -bf N      strategy for bilateral filtering\n");
  printf("              0.no bilateral; 1.independent; 2.composed\n");
  printf("              default value %d\n", PAR_DEFAULT_BILATERAL);
  printf("   -st N      Gaussian standard deviation for temporal dimension\n");
  printf("              default value %f\n", PAR_DEFAULT_SIGMA_T);
  printf("   -s11 N     Gaussian standard deviation for h11\n");
  printf("              default value %f\n", PAR_DEFAULT_SIGMA_H);
  printf("   -s12 N     Gaussian standard deviation for h12\n");
  printf("              default value %f\n", PAR_DEFAULT_SIGMA_H);
  printf("   -s13 N     Gaussian standard deviation for h13\n");
  printf("              default value %f\n", PAR_DEFAULT_SIGMA_H);
  printf("   -s21 N     Gaussian standard deviation for h21\n");
  printf("              default value %f\n", PAR_DEFAULT_SIGMA_H);
  printf("   -s22 N     Gaussian standard deviation for h22\n");
  printf("              default value %f\n", PAR_DEFAULT_SIGMA_H);
  printf("   -s23 N     Gaussian standard deviation for h23\n");
  printf("              default value %f\n", PAR_DEFAULT_SIGMA_H);
  printf("   -s31 N     Gaussian standard deviation for h31\n");
  printf("              default value %f\n", PAR_DEFAULT_SIGMA_H);
  printf("   -s32 N     Gaussian standard deviation for h32\n");
  printf("              default value %f\n", PAR_DEFAULT_SIGMA_H);
  printf("   -b N     type of boundary condition: \n");
  printf("              0.constant; 1.Neumann; 2.Dirichlet\n");
  printf("              default value %d\n", PAR_DEFAULT_BC);
  printf("   -p N     video postprocessing \n");
  printf("              0.no postprocessing; 1.fast crop&zoom; 2.crop&zoom\n");
  printf("              default value %d\n", PAR_DEFAULT_POSTPROCESS);
  printf("   -w name  write transformations to file\n");
  printf("   -l name  load transformations from file\n");
  printf("   -f name  write stabilizing transformations to file\n");
  printf("   -online  switch on estadeo online method\n");
  printf("   -v       switch on verbose mode \n\n\n");
}

/**
 *
 *  Read command line parameters 
 *
 */
int read_parameters(
  int   argc, 
  char  *argv[], 
  char  **video_in,
  char  *video_out,
  char  **out_transform,
  char  **in_transform,
  char  **out_smooth_transform,
  int   &width,
  int   &height,
  int   &nframes,
  int   &nparams,
  int   &online, 
  int   &smooth_strategy,
  int   &bilateral,
  float *sigma,
  int   &boundary_condition, 
  int   &postprocessing,
  int   &verbose
)
{
  if (argc < 5){
    print_help(argv[0]); 
    return 0;
  }
  else{
    int i=1;
    *video_in=argv[i++];
    width=atoi(argv[i++]);
    height=atoi(argv[i++]);
    nframes=atoi(argv[i++]);

    *in_transform=NULL;
    *out_transform=NULL;
    *out_smooth_transform=NULL;
    online=0;
    
    //assign default values to the parameters
    strcpy(video_out,PAR_DEFAULT_OUTVIDEO);
    nparams=PAR_DEFAULT_TRANSFORM;
    smooth_strategy=PAR_DEFAULT_SMOOTHING;
    bilateral=PAR_DEFAULT_BILATERAL;
    sigma[0]=PAR_DEFAULT_SIGMA_T;
    sigma[1]=PAR_DEFAULT_SIGMA_H;
    sigma[2]=PAR_DEFAULT_SIGMA_H;
    sigma[3]=PAR_DEFAULT_SIGMA_H;
    sigma[4]=PAR_DEFAULT_SIGMA_H;
    sigma[5]=PAR_DEFAULT_SIGMA_H;
    sigma[6]=PAR_DEFAULT_SIGMA_H;
    sigma[7]=PAR_DEFAULT_SIGMA_H;
    sigma[8]=PAR_DEFAULT_SIGMA_H;
    boundary_condition=PAR_DEFAULT_BC;
    postprocessing=PAR_DEFAULT_POSTPROCESS;
    verbose=PAR_DEFAULT_VERBOSE;
    
    //read each parameter from the command line
    while(i<argc)
    {
      if(strcmp(argv[i],"-o")==0)
        if(i<argc-1)
          strcpy(video_out,argv[++i]);

      if(strcmp(argv[i],"-t")==0)
        if(i<argc-1)
          nparams=atof(argv[++i]);

      if(strcmp(argv[i],"-m")==0)
        if(i<argc-1)
          smooth_strategy=atoi(argv[++i]);

      if(strcmp(argv[i],"-bf")==0)
        if(i<argc-1)
          bilateral=atoi(argv[++i]);

      if(strcmp(argv[i],"-st")==0)
        if(i<argc-1)
          sigma[0]=atof(argv[++i]);
        
      if(strcmp(argv[i],"-s11")==0)
        if(i<argc-1)
          sigma[1]=atof(argv[++i]);
        
      if(strcmp(argv[i],"-s12")==0)
        if(i<argc-1)
          sigma[2]=atof(argv[++i]);
        
      if(strcmp(argv[i],"-s13")==0)
        if(i<argc-1)
          sigma[3]=atof(argv[++i]);
        
      if(strcmp(argv[i],"-s21")==0)
        if(i<argc-1)
          sigma[4]=atof(argv[++i]);

      if(strcmp(argv[i],"-s22")==0)
        if(i<argc-1)
          sigma[5]=atof(argv[++i]);

      if(strcmp(argv[i],"-s23")==0)
        if(i<argc-1)
          sigma[6]=atof(argv[++i]);

      if(strcmp(argv[i],"-s31")==0)
        if(i<argc-1)
          sigma[7]=atof(argv[++i]);

      if(strcmp(argv[i],"-s32")==0)
        if(i<argc-1)
          sigma[8]=atof(argv[++i]);

      if(strcmp(argv[i],"-b")==0)
        if(i<argc-1)
          boundary_condition=atoi(argv[++i]);

      if(strcmp(argv[i],"-p")==0)
        if(i<argc-1)
          postprocessing=atoi(argv[++i]);

      if(strcmp(argv[i],"-w")==0)
        if(i<argc-1)
          *out_transform=argv[++i];

      if(strcmp(argv[i],"-l")==0)
        if(i<argc-1)
          *in_transform=argv[++i];

      if(strcmp(argv[i],"-f")==0)
        if(i<argc-1)
          *out_smooth_transform=argv[++i];

      if(strcmp(argv[i],"-online")==0)
        online=1;

      if(strcmp(argv[i],"-v")==0)
        verbose=1;
      
      i++;
    }

    //check parameter values
    if(nparams!=2 && nparams!=3 && nparams!=4 && 
       nparams!=6 && nparams!=8) nparams=PAR_DEFAULT_TRANSFORM;
    if(smooth_strategy<0 || smooth_strategy>N_SMOOTH_METHODS)
       smooth_strategy=PAR_DEFAULT_SMOOTHING;
    if(boundary_condition<0 || boundary_condition>2)
       boundary_condition=PAR_DEFAULT_BC;
    if(postprocessing<0 || postprocessing>2)
       postprocessing=PAR_DEFAULT_POSTPROCESS;
    if(bilateral<0 || bilateral>2)
       bilateral=PAR_DEFAULT_BILATERAL;
    if(sigma[0]<0.01)
       sigma[0]=0.01;
    if(sigma[1]<1E-5)
       sigma[1]=1E-5;
    if(sigma[2]<1E-5)
       sigma[2]=1E-5;
    if(sigma[3]<1E-5)
       sigma[3]=1E-5;
    if(sigma[4]<1E-5)
       sigma[4]=1E-5;
    if(sigma[5]<1E-5)
       sigma[5]=1E-5;
    if(sigma[6]<1E-5)
       sigma[6]=1E-5;
    if(sigma[7]<1E-5)
       sigma[7]=1E-5;
    if(sigma[8]<1E-5)
       sigma[8]=1E-5;
  }

  return 1;
}


/**
  *
  *  Function for converting an rgb image to grayscale levels
  * 
**/
void rgb2gray(
  float *rgb,  //input color image
  float *gray, //output grayscale image
  int nx,      //number of pixels
  int ny, 
  int nz
)
{
  int size=nx*ny;
  if(nz>=3)
    #pragma omp parallel for
    for(int i=0;i<size;i++)
      gray[i]=(0.2989*rgb[i*nz]+0.5870*rgb[i*nz+1]+0.1140*rgb[i*nz+2]);
  else
    #pragma omp parallel for
    for(int i=0;i<size;i++)
      gray[i]=rgb[i];
}



/**
 *
 *  Main program:
 *   This program reads the parameters from the console and
 *   then call the video stabilization method
 *
 */
int main (int argc, char *argv[])
{
  //parameters of the method
  char  *video_in, video_out[300];
  char  *out_transform, *in_transform, *out_smooth_transform;
  int   width, height, nchannels=3, nframes;
  int   nparams, online, smooth_strategy; 
  int   boundary_condition, postprocessing, bilateral, verbose;
  float sigma[9];
  
  //read the parameters from the console
  int result=read_parameters(
    argc, argv, &video_in, video_out, &out_transform, &in_transform, 
    &out_smooth_transform, width, height, nframes, nparams, online,
    smooth_strategy, bilateral, sigma, boundary_condition, postprocessing, 
    verbose
  );
  
  if(result)
  {

    if(verbose)
      printf(
        " Input video: '%s'\n Output video: '%s'\n Width: %d, Height: %d "
        " Number of frames: %d\n Transformation: %d, Smoothing: %d, " 
        " bilateral: %d, Sigmas: %f, %f, %f, %f, %f, %f, %f, %f, %f\n"
        " BC: %d Postprocess: %d\n",
        video_in, video_out, width, height, nframes, nparams,
        smooth_strategy, bilateral, sigma[0], sigma[1], sigma[2], sigma[3], 
        sigma[4], sigma[5], sigma[6], sigma[7], sigma[8], boundary_condition, 
        postprocessing
      );
    
    int fsize=width*height;
    int csize=fsize*nchannels;
    int vsize=csize*nframes;
    unsigned char *I=new unsigned char[vsize];
   
    size_t r=read_video(video_in, I, vsize);
    
    if(r<=0)
    {
      fprintf(stderr, "Error: Cannot read the input video '%s'.\n", video_in);
      return EXIT_FAILURE;
    }

    if(verbose) 
      printf(" Size of video in bytes %d\n", (int) r);

    //convert the input video to float and gray levels
    float *Ic=new float[vsize];
    float *Ig=new float[fsize*nframes];
    for(int i=0; i<vsize; i++)
      Ic[i]=(float)I[i];
      
    for(int i=0; i<nframes; i++)
      rgb2gray(&Ic[csize*i], &Ig[fsize*i], width, height, nchannels);
    
    if(verbose)
      printf("\n Starting the stabilization\n");

    //call the method for video stabilization
    if(online)
      estadeo_online(
        Ig, Ic, width, height, nchannels, nframes, nparams, smooth_strategy,
        bilateral, sigma, boundary_condition, postprocessing, 
        out_smooth_transform, verbose
      );
    else
      estadeo(
        Ig, Ic, width, height, nchannels, nframes, nparams, smooth_strategy, 
        bilateral, sigma, boundary_condition, postprocessing, in_transform, 
        out_transform, out_smooth_transform, verbose
      );
    
    //convert the stabilized video to unsigned char
    for(int i=0; i<vsize; i++)
    {
      if(Ic[i]<0)
        I[i]=0;
      else if(Ic[i]>255)
        I[i]=255;
      else I[i]=(unsigned char)Ic[i];
    }
    
    if(verbose)
      printf("\n Writing the output video to '%s'\n", video_out);

    write_video(video_out, I, vsize);
    
    delete []I;
    delete []Ic;
    delete []Ig;
  }

  return EXIT_SUCCESS;
}
