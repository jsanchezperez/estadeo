#ifndef MOTION_SMOOTHING_H
#define MOTION_SMOOTHING_H

#define PURE_COMPOSITION 0
#define COMPOSITIONAL_SMOOTHING 1
#define COMPOSITIONAL_LOCAL_SMOOTHING 2
#define LOCAL_MATRIX_BASED_SMOOTHING 3
#define LOCAL_POINT_BASED_SMOOTHING 4
#define LOCAL_LINEAR_MATRIX_BASED_SMOOTHING 5
#define LOCAL_LINEAR_POINT_BASED_SMOOTHING 6
#define N_SMOOTH_METHODS LOCAL_LINEAR_POINT_BASED_SMOOTHING

//boundary conditions
#define CONSTANT_BC  0
#define NEUMANN_BC   1
#define DIRICHLET_BC 2

//bilateral filters
#define NO_BILATERAL  0
#define BILATERAL_INDEPENDENT 1
#define BILATERAL_COMPOSED 2

void compute_smooth_transforms(
  float *H,       //original matrix transformations
  float *Hp,      //stabilizing transformations
  float *Ho,      //output smooth camera path
  int nparams,    //type of transformation
  int ntransforms //number of frames in the video
);


void motion_smoothing
(
  float *H,          //original matrix transformations
  float *Hp,         //smooth output matrix transformations
  int   nparams,     //type of matrix transformation
  int   ntransforms, //number of frames of the video
  int   bilateral,   //strategies for the bilateral filter
  float *sigma,      //Gaussian standard deviations
  int   type,        //motion smoothing strategy
  int   bc,          //type of boundary condition
  int   verbose      //verbose mode
);


void online_smoothing
(
  float *H,          //original matrix transformations
  float *Hp,         //smooth output matrix transformations
  int   nparams,     //type of matrix transformation
  int   ntransforms, //number of frames of the video
  int   bilateral,   //strategies for the bilateral filter
  float *sigma,      //Gaussian standard deviations
  int   type,        //motion smoothing strategy
  int   bc           //type of boundary condition
);

#endif
