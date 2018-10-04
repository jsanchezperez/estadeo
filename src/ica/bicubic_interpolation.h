// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2015, Javier Sánchez Pérez <jsanchez@ulpgc.es>
// Copyright (C) 2014, Nelson Monzón López <nmonzon@ctim.es>
// All rights reserved.


#ifndef BICUBIC_INTERPOLATION_H
#define BICUBIC_INTERPOLATION_H

#include <vector>

/**
  *
  * Neumann boundary condition test
  *
**/
int neumann_bc (int x, int nx);


/**
  *
  * Compute the bicubic interpolation of a point in an image. 
  * Detects if the point goes outside the image domain
  *
**/
float
bicubic_interpolation(
  float *input,//image to be interpolated
  float uu,    //x component of the vector field
  float vv,    //y component of the vector field
  int nx,       //width of the image
  int ny,       //height of the image
  bool border_out = false //if true, put zeros outside the region
);



/**
  *
  * Bicubic interpolation in one dimension
  *
**/
float
cubic_interpolation(
  float v[4],  //interpolation points
  float x      //point to be interpolated
);



/**
  *
  * Bicubic interpolation in two dimension
  *
**/
float
bicubic_interpolation(
  float p[4][4], //array containing the interpolation points
  float x,       //x position to be interpolated
  float y        //y position to be interpolated
);

/**
  *
  * Compute the bicubic interpolation of an image from a parametric trasform
  *
**/
void bicubic_interpolation(
  float *input,        //image to be warped
  std::vector<int> &x, //selected points
  float *output,       //warped output image with bicubic interpolation
  float *params,       //x component of the vector field
  int nparams,          //number of parameters of the transform
  int nx,               //width of the image
  int ny,               //height of the image
  bool border_out=true  //if true, put zeros outside the region
);



/**
  *
  * Compute the bicubic interpolation of an image from a parametric trasform
  *
**/
void bicubic_interpolation(
  float *input,        //image to be warped
  float *output,       //warped output image with bicubic interpolation
  float *params,       //x component of the vector field
  int nparams,          //number of parameters of the transform
  int nx,               //width of the image
  int ny,               //height of the image
  bool border_out=true  //if true, put zeros outside the region
);




#endif
