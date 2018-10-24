// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2017-2018, Javier Sánchez Pérez <jsanchez@ulpgc.es>
// All rights reserved.


#include "motion_smoothing.h"
#include "gaussian_conv_dct.h"
#include "transformation.h"
#include "matrix.h"

#include <stdio.h>
#include <math.h>


//Gaussian convolution
void online_global_gaussian
(
  float *H,          //original matrix transformations
  float *Hs,         //smooth output matrix transformations
  int   i,           //frame number
  int   nparams,     //type of matrix transformation
  int   ntransforms, //number of frames of the video  
  float sigma        //Gaussian standard deviation

)
{
  if(sigma>=3*ntransforms)
    sigma=ntransforms/3;
  
  int radius=3*sigma;

  if(radius>=ntransforms) 
    radius=ntransforms-1;
  
  //Gaussian convolution in  each parameter separately
  for(int p=0;p<nparams;p++)
  {
    double average=0.0;
    double sum=0.0;
    
    for(int j=i-radius;j<=i+radius;j++)
    {
      double value=0;
      
      //Neumann boundary conditions
      if(j<0)
        value=H[-j*nparams+p];      
      else if(j>=ntransforms)
        value=H[(2*ntransforms-1-j)*nparams+p];
      else 
        value=H[j*nparams+p];
      
      //increase accumulator
      double norm=0.5*(j-i)*(j-i)/(sigma*sigma);
      double gauss=exp(-norm);
      average+=gauss*value;
      sum+=gauss;
    }
    Hs[p]=(float) (average/sum);
  }
}


//Gaussian convolution for local methods
void online_local_gaussian
(
  float *H,          //original matrix transformations
  float *Hs,         //smooth output matrix transformations
  int   i,           //frame number
  int   nparams,     //type of matrix transformation
  int   ntransforms, //number of frames of the video  
  float sigma,       //Gaussian standard deviation
  float *H_1         //inverse transforms
)
{
  if(sigma>=3*ntransforms)
    sigma=ntransforms/3;
  
  int radius=3*sigma;

  if(radius>=ntransforms) 
    radius=ntransforms-1;
  
  //Gaussian convolution in each parameter separately
  for(int p=0;p<nparams;p++)
  { 
    float average=0.0;
    float sum=0.0;
    
    for(int j=i-radius;j<=i+radius;j++)
    {
      float value=0;
      
      //test boundary conditions
      if(j<0)
        value=H_1[(-j)*nparams+p];
      else if(j>=ntransforms) 
        value=H_1[(2*ntransforms-1-j)*nparams+p];
      else 
        value=H[j*nparams+p];
      
      //increase accumulator
      float norm=0.5*(j-i)*(j-i)/(sigma*sigma);
      float gauss=exp(-norm);
      average+=gauss*value;
      sum+=gauss;
    }
    Hs[p]=average/sum;
  }
}


//local matrix based smoothing approach
void online_local_matrix_based_smoothing
(
  float *H,          //original matrix transformations
  float *Hp,         //smooth output matrix transformations
  int   nparams,     //type of matrix transformation
  int   ntransforms, //number of frames of the video
  float sigma        //Gaussian standard deviation
)
{
  float *H_1=new float[ntransforms*nparams];
  float *Hc=new float[ntransforms*nparams];
  float *Hs=new float[ntransforms*nparams];
  
  float sigma_0=sigma[0];
  if(sigma_0>=3*ntransforms)
    sigma_0=ntransforms/3;
  
  int radius=3*sigma_0;

  if(radius>=ntransforms) 
    radius=ntransforms-1;

  //compute inverse transformations 
  for(int i=0;i<ntransforms;i++) 
    inverse_transform(&(H[i*nparams]), &(H_1[i*nparams]), nparams);

  //recompute the stabilization for past frames
  for(int i=ntransforms-radius;i<ntransforms;i++)
  {
    int t1=(i-radius>0)?i-radius:0;
    int t2=((i+radius)<ntransforms)?(i+radius):ntransforms-1;

    //compute backward transformations
    if(i>0)
    {
      for(int j=0;j<nparams;j++) 
         Hc[(i-1)*nparams+j]=H_1[i*nparams+j];
      for(int j=i-2;j>=t1;j--)
         compose_transform(
            &(H_1[(j+1)*nparams]), &(Hc[(j+1)*nparams]), 
            &(Hc[j*nparams]), nparams
         );
    }

    //introduce the identity matrix in the middle
    for(int j=0;j<nparams;j++) Hc[i*nparams+j]=0;

    //compute forward transformations
    if(i<ntransforms-1)
    {
      for(int j=0;j<nparams;j++) 
         Hc[(i+1)*nparams+j]=H[(i+1)*nparams+j];
      for(int j=i+2;j<=t2;j++)
         compose_transform(
            &(H[j*nparams]), &(Hc[(j-1)*nparams]), &(Hc[j*nparams]), nparams
         );
    }

    //smooth transforms with a discrete Gaussian kernel
    global_gaussian(
      H, Hc, &(Hs[i*nparams]), i, nparams, ntransforms, sigma
    );

    //compute inverse transformations 
    inverse_transform(&(Hs[i*nparams]), &(Hp[i*nparams]), nparams);
  }
  
  delete []H_1;  
  delete []Hc;  
  delete []Hs;
}




//Gaussian convolution with a set of points
float online_point_gaussian(
  float *x,          //set of points
  int   i,           //frame number
  int   ntransforms, //number of transforms
  float sigma        //Gaussian standard deviation
)
{
  float average=0.0;
  float sum=0.0;

  if(sigma>=3*ntransforms)
    sigma=ntransforms/3;
  
  int radius=3*sigma;

  if(radius>=ntransforms) 
    radius=ntransforms-1;
  
  //Gaussian convolution
  for(int j=i-radius;j<=i+radius;j++)
  {
    float value=0;
    
    //test boundary conditions
    if(j<0)
      value=x[-j];
    else if(j>=ntransforms) 
      value=x[2*ntransforms-1-j];
    else 
      value=x[j];

    float dx=j-i;
    float norm=0.5*dx*dx/(sigma*sigma);
    float gauss=exp(-norm);
    average+=gauss*value;
    sum+=gauss;
  }
  return average/sum;
}


//Matrix DCT Gaussian convolution
void online_matrix_gaussian_dct
(
  float *H,      //original matrix transformations
  float *Hs,     //smooth output matrix transformations
  int   nparams, //type of matrix transformation
  int   N,       //number of frames of the video  
  float sigma    //Gaussian standard deviation
)
{  
  num *dest=new num[3*N-2];
  num *src=new num[3*N-2];
  dct_coeffs c;
  
  //convolution in each matrix position
  for(int p=0; p<nparams; p++)
  {
    //copy the original image
    for(int i=N-1; i<2*N-1; i++)
      src[i]=H[(i-N+1)*nparams+p];
    
    //Neumann boundary conditions
    for(int i=0; i<N-1; i++)
      src[i]=H[(N-i-2)*nparams+p];
    for(int i=2*N-1; i<3*N-2; i++)
      src[i]=H[(3*N-i-2)*nparams+p];

    //apply DCT Gaussian convolution
    if (!(dct_precomp(&c, dest, src, 3*N-2, 1, sigma)))
      printf("Error in Gaussian convolution with DCT.");
    else {
      dct_gaussian_conv(c);
      dct_free(&c);
    }

    //copy the signal in the domain
    for(int i=0; i<N; i++)
      Hs[i*nparams+p]=dest[N-1+i];
  }

  delete []src;
  delete []dest;
}


//local linear matrix-based smoothing
void online_local_linear_matrix_based_smoothing
(
  float *H,          //original matrix transformations
  float *Hp,         //smooth output matrix transformations
  int   nparams,     //type of matrix transformation
  int   ntransforms, //number of frames of the video  
  float sigma        //Gaussian standard deviations
)
{
  float *HH=new float[ntransforms*9];
  float *Hi=new float[ntransforms*9];
  float *Hs=new float[ntransforms*9];

  //convert from params to matrices
  for(int i=1;i<ntransforms;i++) 
    params2matrix(&(H[i*nparams]), &(HH[i*9]), nparams);

  //identity matrix in the first position
  Hi[1]=Hi[2]=Hi[3]=0;
  Hi[5]=Hi[6]=Hi[7]=0;
  Hi[0]=Hi[4]=Hi[8]=1;

  //compute the virtual matrix trajectories
  for(int i=1;i<ntransforms;i++) 
  {
    for(int j=0;j<9;j++)
      //accumulate the homographies
      Hi[i*9+j]=Hi[(i-1)*9+j]+HH[i*9+j];
    
    //subtract the identity matrix
    Hi[i*9]-=1;Hi[i*9+4]-=1;Hi[i*9+8]-=1;
  }

  //convolve the virtual trajectories with a Gaussian kernel
  online_matrix_gaussian_dct(Hi, Hs, 9, ntransforms, sigma);

  for(int i=0;i<ntransforms;i++) 
  {
    //compute the correction matrix
    for(int j=0;j<9;j++)
      Hs[i*9+j]-=Hi[i*9+j];

    //add the identity matrix
    Hs[i*9]+=1;Hs[i*9+4]+=1;Hs[i*9+8]+=1;

    //convert homographies to params
    matrix2params(&(Hs[i*9]),&(Hp[i*nparams]),nparams);

    //compute its inverse 
    inverse_transform(&(Hp[i*nparams]), &(Hp[i*nparams]), nparams);        
  }
  
  delete []HH;
  delete []Hi;
  delete []Hs;
}



//Point DCT Gaussian convolution
void online_point_gaussian_dct
(
  float *x,    //input set of points
  float *xs,   //output set of smoothed points 
  int   N,     //number of points
  float sigma, //Gaussian standard deviation
  int   bc     //type of boundary condition
)
{  
  num *dest=new num[3*N-2];
  num *src=new num[3*N-2];
  dct_coeffs c;

  //copy the original signal
  for(int i=N-1; i<2*N-1; i++)
    src[i]=x[i-N+1];
  
  //boundary conditions
  switch(bc){
    case CONSTANT_BC:
      for(int i=0; i<N-1; i++)
        src[i]=x[0];
      for(int i=2*N-1; i<3*N-2; i++)
        src[i]=x[N-1];
      break;
    case NEUMANN_BC:
      for(int i=0; i<N-1; i++)
        src[i]=x[N-i-2];
      for(int i=2*N-1; i<3*N-2; i++)
        src[i]=x[3*N-i-2];
      break;
    case DIRICHLET_BC: 
      for(int i=0; i<N-1; i++)
        src[i]=2*x[0]-x[N-i-2];
      for(int i=2*N-1; i<3*N-2; i++)
        src[i]=2*x[N-1]-x[3*N-i-2];
      break;
  }

  //apply DCT Gaussian convolution
  if (!(dct_precomp(&c, dest, src, 3*N-2, 1, sigma)))
    printf("Error in Gaussian convolution with DCT.");
  else {
    dct_gaussian_conv(c);
    dct_free(&c);
  }

  //copy the signal in the domain
  for(int i=0; i<N; i++)
    xs[i]=dest[N-1+i];
  
  delete []src;
  delete []dest;
}


//local linear point based smoothing approach
void online_local_linear_point_based_smoothing
(
  float *H,          //original matrix transformations
  float *Hp,         //smooth output matrix transformations
  int   nparams,     //type of matrix transformation
  int   ntransforms, //number of frames of the video  
  float sigma        //Gaussian standard deviations
)
{ 
  float *x0=new  float[ntransforms];
  float *x1=new  float[ntransforms];
  float *x2=new  float[ntransforms];
  float *x3=new  float[ntransforms];
  float *y0=new  float[ntransforms];
  float *y1=new  float[ntransforms];
  float *y2=new  float[ntransforms];
  float *y3=new  float[ntransforms];
  float *x0s=new float[ntransforms];
  float *x1s=new float[ntransforms];
  float *x2s=new float[ntransforms];
  float *x3s=new float[ntransforms];
  float *y0s=new float[ntransforms];
  float *y1s=new float[ntransforms];
  float *y2s=new float[ntransforms];
  float *y3s=new float[ntransforms];
  float *HH=new  float[ntransforms*9];

  //convert from params to matrices
  for(int i=0;i<ntransforms;i++) 
    params2matrix(&(H[i*nparams]), &(HH[i*9]), nparams);

  //choose four fixed points for all frames
  float xp[4]={0, 0, 500, 500};
  float yp[4]={0, 500, 0, 500};

  //tracking a set of points to be smoothed
  x0[0]=xp[0]; y0[0]=yp[0];
  x1[0]=xp[1]; y1[0]=yp[1];
  x2[0]=xp[2]; y2[0]=yp[2];
  x3[0]=xp[3]; y3[0]=yp[3];
  
  //compute the virtual trajectories
  for(int i=1;i<ntransforms;i++) 
  { 
    float dx, dy;
    Hx(&(HH[i*9]),xp[0],yp[0],dx,dy);
    x0[i]=x0[i-1]+(dx-xp[0]);
    y0[i]=y0[i-1]+(dy-yp[0]);

    Hx(&(HH[i*9]),xp[1],yp[1],dx,dy);
    x1[i]=x1[i-1]+(dx-xp[1]);
    y1[i]=y1[i-1]+(dy-yp[1]);

    Hx(&(HH[i*9]),xp[2],yp[2],dx,dy);
    x2[i]=x2[i-1]+(dx-xp[2]);
    y2[i]=y2[i-1]+(dy-yp[2]);

    Hx(&(HH[i*9]),xp[3],yp[3],dx,dy);
    x3[i]=x3[i-1]+(dx-xp[3]);
    y3[i]=y3[i-1]+(dy-yp[3]);
  }
  
  //DCT Gaussian convolution of each virtual trajectory
  online_point_gaussian_dct(x0, x0s, ntransforms, sigma);
  online_point_gaussian_dct(x1, x1s, ntransforms, sigma);
  online_point_gaussian_dct(x2, x2s, ntransforms, sigma);
  online_point_gaussian_dct(x3, x3s, ntransforms, sigma);
  online_point_gaussian_dct(y0, y0s, ntransforms, sigma);
  online_point_gaussian_dct(y1, y1s, ntransforms, sigma);
  online_point_gaussian_dct(y2, y2s, ntransforms, sigma);
  online_point_gaussian_dct(y3, y3s, ntransforms, sigma);
  
  for(int i=0;i<ntransforms;i++) 
  {
    x0s[i]+=xp[0]-x0[i];
    x1s[i]+=xp[1]-x1[i];
    x2s[i]+=xp[2]-x2[i];
    x3s[i]+=xp[3]-x3[i];
    y0s[i]+=yp[0]-y0[i];
    y1s[i]+=yp[1]-y1[i];
    y2s[i]+=yp[2]-y2[i];
    y3s[i]+=yp[3]-y3[i];
    
    //calculate the smoothed homography    
    float tmp[9];
    compute_H(
      x0s[i],x1s[i],x2s[i],x3s[i],
      y0s[i],y1s[i],y2s[i],y3s[i],      
      xp[0],xp[1],xp[2],xp[3],      
      yp[0],yp[1],yp[2],yp[3],     
      tmp
    );

    float Hout[9];
    for(int j=0; j<9; j++) Hout[j]=tmp[j];
    
    //convert homographies to params
    matrix2params(Hout,&(Hp[i*nparams]),nparams);
  }

  delete []x0s;
  delete []x1s;
  delete []x2s;
  delete []x3s;
  delete []y0s;
  delete []y1s;
  delete []y2s;
  delete []y3s;
  delete []x0;
  delete []x1;
  delete []x2;
  delete []x3;
  delete []y0;
  delete []y1;
  delete []y2;
  delete []y3;
  delete []HH;
}


/**
  *
  * Main function for motion_smoothing
  * 
  *
**/
void online_smoothing
(
  float *H,          //original matrix transformations
  float *Hp,         //smooth output matrix transformations
  int   nparams,     //type of matrix transformation
  int   ntransforms, //number of frames of the video
  float sigma,       //Gaussian standard deviations
  int   type         //motion smoothing strategy
)
{
  switch(type) 
  {
    default: case LOCAL_MATRIX_BASED_SMOOTHING:
      online_local_matrix_based_smoothing(
        H, Hp, nparams, ntransforms, bilateral, sigma, bc
      );
      break;
    case LOCAL_LINEAR_MATRIX_BASED_SMOOTHING:
      online_local_linear_matrix_based_smoothing(
        H, Hp, nparams, ntransforms, sigma, bc
      );
      break;
    case LOCAL_LINEAR_POINT_BASED_SMOOTHING:
      online_local_linear_point_based_smoothing(
        H, Hp, nparams, ntransforms, sigma, bc
      );
      break;
  }
}

