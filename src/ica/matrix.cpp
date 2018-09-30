// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2015-2018, Javier Sánchez Pérez <jsanchez@dis.ulpgc.es>
// All rights reserved.

#include "matrix.h"

#include <math.h>

//Multiplication of a square matrix and a vector
void Axb(float *A, float *b, float *p, int n)
{
  for(int i=0; i<n; i++)
  {
    float sum=0;
    for(int j=0; j<n; j++)
      sum+=A[i*n+j]*b[j];
    
    p[i]=sum;
  }   
}




//Multiplication of a square matrix and a vector
void Axb(double *A, double *b, double *p, int n)
{
  for(int i=0; i<n; i++)
  {
    float sum=0;
    for(int j=0; j<n; j++)
      sum+=A[i*n+j]*b[j];
    
    p[i]=sum;
  }   
}

//Function to compute the inverse of a matrix
//through Gaussian elimination
int inverse(
  float *A,   //input matrix
  float *A_1, //output matrix
  int N        //matrix dimension
) 
{
  double *PASO=new double[2*N*N];

  double max,paso,mul;
  int i,j,i_max,k;

  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      PASO[i*2*N+j]=A[i*N+j];
      PASO[i*2*N+j+N]=0.;
    }
  }    
  for(i=0;i<N;i++)
      PASO[i*2*N+i+N]=1.;      
      
  for(i=0;i<N;i++){
    max=fabs(PASO[i*2*N+i]);
    i_max=i;
    for(j=i;j<N;j++){
       if(fabs(PASO[j*2*N+i])>max){
         i_max=j; max=fabs(PASO[j*2*N+i]);
       } 
    }

    if(max<10e-30){ 
      delete []PASO;
      return -1;
    }
    if(i_max>i){
      for(k=0;k<2*N;k++){
        paso=PASO[i*2*N+k];
        PASO[i*2*N+k]=PASO[i_max*2*N+k];
        PASO[i_max*2*N+k]=paso;
      }
    } 

    for(j=i+1;j<N;j++){
      mul=-PASO[j*2*N+i]/PASO[i*2*N+i];
      for(k=i;k<2*N;k++) PASO[j*2*N+k]+=mul*PASO[i*2*N+k];                
    }
  }
  
  if(fabs(PASO[(N-1)*2*N+N-1])<10e-30){ 
      delete []PASO;
      return -1;
  }
      
  for(i=N-1;i>0;i--){
    for(j=i-1;j>=0;j--){
      mul=-PASO[j*2*N+i]/PASO[i*2*N+i];
      for(k=i;k<2*N;k++) PASO[j*2*N+k]+=mul*PASO[i*2*N+k];     
    }
  }  
  for(i=0;i<N;i++)
    for(j=N;j<2*N;j++)
      PASO[i*2*N+j]/=PASO[i*2*N+i];  
    
  for(i=0;i<N;i++)
    for(j=0;j<N;j++)
      A_1[i*N+j]=PASO[i*2*N+j+N];

  delete []PASO;
  
  return 0;   
}


//Function to compute the inverse of a matrix
//through Gaussian elimination
int inverse(
  double *A,   //input matrix
  double *A_1, //output matrix
  int N        //matrix dimension
) 
{
  double *PASO=new double[2*N*N];

  double max,paso,mul;
  int i,j,i_max,k;

  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      PASO[i*2*N+j]=A[i*N+j];
      PASO[i*2*N+j+N]=0.;
    }
  }    
  for(i=0;i<N;i++)
      PASO[i*2*N+i+N]=1.;      
      
  for(i=0;i<N;i++){
    max=fabs(PASO[i*2*N+i]);
    i_max=i;
    for(j=i;j<N;j++){
       if(fabs(PASO[j*2*N+i])>max){
         i_max=j; max=fabs(PASO[j*2*N+i]);
       } 
    }

    if(max<10e-30){ 
      delete []PASO;
      return -1;
    }
    if(i_max>i){
      for(k=0;k<2*N;k++){
        paso=PASO[i*2*N+k];
        PASO[i*2*N+k]=PASO[i_max*2*N+k];
        PASO[i_max*2*N+k]=paso;
      }
    } 

    for(j=i+1;j<N;j++){
      mul=-PASO[j*2*N+i]/PASO[i*2*N+i];
      for(k=i;k<2*N;k++) PASO[j*2*N+k]+=mul*PASO[i*2*N+k];                
    }
  }
  
  if(fabs(PASO[(N-1)*2*N+N-1])<10e-30){ 
      delete []PASO;
      return -1;
  }
      
  for(i=N-1;i>0;i--){
    for(j=i-1;j>=0;j--){
      mul=-PASO[j*2*N+i]/PASO[i*2*N+i];
      for(k=i;k<2*N;k++) PASO[j*2*N+k]+=mul*PASO[i*2*N+k];     
    }
  }  
  for(i=0;i<N;i++)
    for(j=N;j<2*N;j++)
      PASO[i*2*N+j]/=PASO[i*2*N+i];  
    
  for(i=0;i<N;i++)
    for(j=0;j<N;j++)
      A_1[i*N+j]=PASO[i*2*N+j+N];

  delete []PASO;
  
  return 0;   
}


//Multiplies a 3x3 homography and a vector
void Hx(float *H, float x, float y, float &xp, float &yp)
{
  float den=H[6]*x+H[7]*y+H[8];
  if(den*den>1E-10)
  {
    xp=(H[0]*x+H[1]*y+H[2])/den;
    yp=(H[3]*x+H[4]*y+H[5])/den;
  }
  else {
    xp=x;yp=y;
  }
}

//Multiplies two 3x3 matrices
void HxH(float *H1, float *H2, float *H3)
{
  H3[0]=H1[0]*H2[0]+H1[1]*H2[3]+H1[2]*H2[6];
  H3[1]=H1[0]*H2[1]+H1[1]*H2[4]+H1[2]*H2[7];
  H3[2]=H1[0]*H2[2]+H1[1]*H2[5]+H1[2]*H2[8];
  H3[3]=H1[3]*H2[0]+H1[4]*H2[3]+H1[5]*H2[6];
  H3[4]=H1[3]*H2[1]+H1[4]*H2[4]+H1[5]*H2[7];
  H3[5]=H1[3]*H2[2]+H1[4]*H2[5]+H1[5]*H2[8];
  H3[6]=H1[6]*H2[0]+H1[7]*H2[3]+H1[8]*H2[6];
  H3[7]=H1[6]*H2[1]+H1[7]*H2[4]+H1[8]*H2[7];
  H3[8]=H1[6]*H2[2]+H1[7]*H2[5]+H1[8]*H2[8];
}


//Computes a homography from four points
void compute_H
(
  double x1, double x2, double x3, double x4,          
  double y1, double y2, double y3, double y4,
  double x1p, double x2p, double x3p, double x4p,          
  double y1p, double y2p, double y3p, double y4p,
  double *H
)
{
  double A[64], A_1[64], b[8];
  
  //compose the independent vector
  b[0]=x1p; b[1]=y1p;
  b[2]=x2p; b[3]=y2p;
  b[4]=x3p; b[5]=y3p;
  b[6]=x4p; b[7]=y4p;

  //compose the system matrix
  int i=0; 
  int j=0;
  A[i++]=x1;A[i++]=y1;A[i++]=1;A[i++]=0;A[i++]=0;
  A[i++]=0;A[i++]=-b[j]*x1;A[i++]=-b[j++]*y1;
  A[i++]=0;A[i++]=0;A[i++]=0;A[i++]=x1;A[i++]=y1;
  A[i++]=1;A[i++]=-b[j]*x1;A[i++]=-b[j++]*y1;
  A[i++]=x2;A[i++]=y2;A[i++]=1;A[i++]=0;A[i++]=0;
  A[i++]=0;A[i++]=-b[j]*x2;A[i++]=-b[j++]*y2;
  A[i++]=0;A[i++]=0;A[i++]=0;A[i++]=x2;A[i++]=y2;
  A[i++]=1;A[i++]=-b[j]*x2;A[i++]=-b[j++]*y2;
  A[i++]=x3;A[i++]=y3;A[i++]=1;A[i++]=0;A[i++]=0;
  A[i++]=0;A[i++]=-b[j]*x3;A[i++]=-b[j++]*y3;
  A[i++]=0;A[i++]=0;A[i++]=0;A[i++]=x3;A[i++]=y3;
  A[i++]=1;A[i++]=-b[j]*x3;A[i++]=-b[j++]*y3;
  A[i++]=x4;A[i++]=y4;A[i++]=1;A[i++]=0;A[i++]=0;
  A[i++]=0;A[i++]=-b[j]*x4;A[i++]=-b[j++]*y4;
  A[i++]=0;A[i++]=0;A[i++]=0;A[i++]=x4;A[i++]=y4;
  A[i++]=1;A[i++]=-b[j]*x4;A[i++]=-b[j++]*y4;

  //solve
  inverse(A, A_1, 8);
  Axb(A_1,b,H,8);
  H[8]=1; 
}





