#include <iostream>
#include <complex>

using namespace std;

#include "Matrix.h"


Matrix::Matrix(const Matrix &m)
{

  for(int i=0; i<4; i++){
    for(int j=0; j<4; j++){
      M[i][j] = m.M[i][j];
    } 
  }  

}


Matrix::Matrix()
{

}


Matrix::Matrix(complex<double> A[4][4])
{
  
  
  for(int i=0; i<4; i++){
    for(int j=0; j<4; j++){
      M[i][j] = A[i][j];
    }
  }
  

}


Matrix::~Matrix()
{


}


const Matrix &Matrix::operator=( const Matrix &right)
{

  if ( &right != this)
    {

      for(int i=0; i<4; i++){
	for(int j=0; j<4; j++){
	  M[i][j] = right.M[i][j];
    	}
      }

    }

  return *this;

}


Matrix ConjMatrix(Matrix m)
{
  
  Matrix temp;

  for(int i=0; i<4; i++){
    for(int j=0; j<4; j++){ 
      temp.M[i][j] = conj(m.M[j][i]); 
    }
  }

  return temp; 

}

Matrix GammaMatrixGamma(Matrix m)
{

  Matrix temp;

  temp.M[0][0] = m.M[0][0];
  temp.M[0][1] = m.M[0][1];
  temp.M[0][2] = -m.M[0][2];
  temp.M[0][3] = -m.M[0][3];

  temp.M[1][0] = m.M[1][0];
  temp.M[1][1] = m.M[1][1];
  temp.M[1][2] = -m.M[1][2];
  temp.M[1][3] = -m.M[1][3];

  temp.M[2][0] = -m.M[2][0];
  temp.M[2][1] = -m.M[2][1];
  temp.M[2][2] = m.M[2][2];
  temp.M[2][3] = m.M[2][3];

  temp.M[3][0] = -m.M[3][0];
  temp.M[3][1] = -m.M[3][1];
  temp.M[3][2] = m.M[3][2];
  temp.M[3][3] = m.M[3][3];

  return temp;

}


complex<double> Trace(Matrix A){

  complex<double> temp;

  for(int i=0; i<4; i++){
    
    temp += A.M[i][i]; 
    
  }

  return temp;

}

void Trace_product(const Matrix &M1, const Matrix &M2, complex<double> &result)
{
	//We use this function to remove some matrix multiplication, where 4*4 loops would be done to do the product, and then an additional 4 to compute the trace of the result
 	for (int i = 0 ; i < 4 ; i++)
	{
		for (int j=0 ; j < 4 ; j++)
		{
			result += M1.M[i][j]*M2.M[j][i];
		}
	} 
}


Matrix operator*(double op1, Matrix op2)
{

  Matrix temp;
  
  for(int i=0; i<4; i++){
    for(int j=0; j<4; j++){ 
      temp.M[i][j] = op1*op2.M[i][j]; 
    }
  }

  return temp;

}


Matrix operator*(complex<double> op1, Matrix op2)
{

  Matrix temp;
  
  for(int i=0; i<4; i++){
    for(int j=0; j<4; j++){ 
      temp.M[i][j] = op1*op2.M[i][j]; 
    }
  }

  return temp;

}


Matrix operator*(Matrix op1, Matrix op2)
{

  Matrix temp;
  
  for(int i=0; i<4; i++){
    for(int j=0; j<4; j++){ 
      for(int k=0; k<4; k++){
      temp.M[i][j] += op1.M[i][k]*op2.M[k][j];
      } 
    }
  }

  return temp;

}


Matrix operator+(Matrix op1, Matrix op2)
{

  Matrix temp;
  
  for(int i=0; i<4; i++){
    for(int j=0; j<4; j++){
      temp.M[i][j] = op1.M[i][j] + op2.M[i][j]; 
    }
  }

  return temp;

}


Matrix operator-(Matrix op1, Matrix op2)
{

  Matrix temp;
  
  for(int i=0; i<4; i++){
    for(int j=0; j<4; j++){
      temp.M[i][j] = op1.M[i][j] - op2.M[i][j]; 
    }
  }

  return temp;

}


ostream &operator<<(ostream &stream, Matrix m)
{

  for(int i=0; i<4; i++){
    stream << "|" << " " ;
    for(int j=0; j<4; j++){
      stream << m.M[i][j] << " " ;
    }
    stream << "|" << endl;
  }  

  return stream;

}

