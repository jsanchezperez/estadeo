// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2018, Thibaud Briand <thibaud.briand@enpc.fr>
// Copyright (C) 2015, Javier Sánchez Pérez <jsanchez@dis.ulpgc.es>
// All rights reserved.

#ifndef MATRIX_H
#define MATRIX_H

//Multiplication of a square matrix and a vector
void Axb(float *A, float *b, float *p, int n);

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

//Multiplication with a scalar
//B should be initialized to zero outside
void sA(float s, float *A, float *B, int m);

//Function to compute the inverse of a matrix
//through Gaussian elimination
int inverse(float *A, float *A_1, int N = 3);

//Multiplies a 3x3 homography and a vector
void Hx(float *H, float x, float y, float &xp, float &yp);

//Multiplies two 3x3 matrices
void HxH(float *H1, float *H2, float *H3);

//Computes a homography from four points
void compute_H
(
  float x1, float x2, float x3, float x4,          
  float y1, float y2, float y3, float y4,
  float x1p, float x2p, float x3p, float x4p,          
  float y1p, float y2p, float y3p, float y4p,
  float *H
);



#endif
