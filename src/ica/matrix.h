// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2015, Javier Sánchez Pérez <jsanchez@dis.ulpgc.es>
// All rights reserved.

#ifndef MATRIX_H
#define MATRIX_H

//Multiplication of a square matrix and a vector
void Axb(float *A, float *b, float *p, int n);

//Multiplication of a square matrix and a vector
void Axb(double *A, double *b, double *p, int n);

//Multiplication of the transpose of a matrix and a vector
//p should be initialized to zero outside
void Atb(float *A, float *b, float *p, int n, int m);

//Multiplication of the transpose of a matrix, a vector and a scalar
//p should be initialized to zero outside
void sAtb(float s, float *A, float *b, float *p, int n, int m);

//Multiplication of the transpose of a matrix and itself
//B should be initialized to zero outside
void AtA(float *A, float *B, int n, int m);

//Multiplication of the transpose of a matrix and itself with a scalar
//B should be initialized to zero outside
void sAtA(float s, float *A, float *B, int n, int m);

//Function to compute the inverse of a matrix
//through Gaussian elimination
int inverse(float *A, float *A_1, int N = 3);

//Function to compute the inverse of a matrix
//through Gaussian elimination
int inverse(double *A, double *A_1, int N = 3);


//Multiplies a 3x3 homography and a vector
void Hx(float *H, float x, float y, float &xp, float &yp);

//Multiplies two 3x3 matrices
void HxH(float *H1, float *H2, float *H3);

//Computes a homography from four points
void compute_H
(
  double x1, double x2, double x3, double x4,          
  double y1, double y2, double y3, double y4,
  double x1p, double x2p, double x3p, double x4p,          
  double y1p, double y2p, double y3p, double y4p,
  double *H
);


#endif
