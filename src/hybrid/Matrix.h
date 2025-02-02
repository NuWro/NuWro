#ifndef MATRIX_H     
#define MATRIX_H

#include <complex>
using std::complex;
#include <iostream>
using std::ostream;
/********************************************************************

The Matrix class stores (4, 4) complex matrices

PRIVATE
the viewpoint we take is that of a matrix M represented by a pointer
to pointers (see web-link). a "first level" pointer is allocated; it
points to an array of #1 pointers, one for each row. next, we allocate
each row-pointer to point to an array of #2 elements. in general, #1
stands for nrows and #2 for ncolumns.

PUBLIC
a copy-constructor is defined, next to a constructor with and w/o 
parameters.      
                                                                    
*********************************************************************/
//A: Matrix has been changed from persistent memory assignment to stack memory and fixed size of 4x4. 
//Much faster but with some downsides when it comes to user-friendliness.

class Matrix{
  friend Matrix ConjMatrix( Matrix );
  friend Matrix GammaMatrixGamma( Matrix );
  friend complex<double> Trace( Matrix );

  friend void Trace_product(const Matrix &M1, const Matrix &M2, complex<double> &result);
  friend Matrix operator*( double, Matrix );
  friend Matrix operator*( complex<double>, Matrix );
  friend Matrix operator*( Matrix, Matrix );
  friend Matrix operator+( Matrix, Matrix );
  friend Matrix operator-( Matrix, Matrix );
  friend ostream &operator<<( ostream &, Matrix );


 public:
  Matrix( const Matrix & );  // copy constructor
  Matrix( complex<double> A[4][4] );     
  Matrix();
  ~Matrix();

  const Matrix &operator=( const Matrix & );

//  private:
  complex<double> M[4][4]={};
};



void Construct_GammaMatrices( Matrix &, Matrix [] );

#endif
