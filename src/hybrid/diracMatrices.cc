// #include <iostream>
// #include <complex>
// 
// using namespace std;

#include "diracMatrices.h"
#include "Constants.h"


complex<double> id[4][4] = { {complex<double> (1, 0), complex<double> (0, 0), complex<double> (0, 0), complex<double> (0, 0)},
			    {complex<double> (0, 0), complex<double> (1, 0), complex<double> (0, 0), complex<double> (0, 0)},
			    {complex<double> (0, 0), complex<double> (0, 0), complex<double> (1, 0), complex<double> (0, 0)},
			    {complex<double> (0, 0), complex<double> (0, 0), complex<double> (0, 0), complex<double> (1, 0)} };

complex<double> h0[4][4] = { {complex<double> (1, 0), complex<double> (0, 0), complex<double> (0, 0), complex<double> (0, 0)},
			    {complex<double> (0, 0), complex<double> (1, 0), complex<double> (0, 0), complex<double> (0, 0)},
			    {complex<double> (0, 0), complex<double> (0, 0), complex<double> (-1, 0), complex<double> (0, 0)},
			    {complex<double> (0, 0), complex<double> (0, 0), complex<double> (0, 0), complex<double> (-1, 0)} };

complex<double> h1[4][4] = { {complex<double> (0, 0), complex<double> (0, 0), complex<double> (0, 0), complex<double> (1, 0)},
			    {complex<double> (0, 0), complex<double> (0, 0), complex<double> (1, 0), complex<double> (0, 0)},
			    {complex<double> (0, 0), complex<double> (-1, 0), complex<double> (0, 0), complex<double> (0, 0)},
			    {complex<double> (-1, 0), complex<double> (0, 0), complex<double> (0, 0), complex<double> (0, 0)} };

complex<double> h2[4][4] = { {complex<double> (0, 0), complex<double> (0, 0), complex<double> (0, 0), complex<double> (0, -1)},
			    {complex<double> (0, 0), complex<double> (0, 0), complex<double> (0, 1), complex<double> (0, 0)},
			    {complex<double> (0, 0), complex<double> (0, 1), complex<double> (0, 0), complex<double> (0, 0)},
			    {complex<double> (0, -1), complex<double> (0, 0), complex<double> (0, 0), complex<double> (0, 0)} };

complex<double> h3[4][4] = { {complex<double> (0, 0), complex<double> (0, 0), complex<double> (1, 0), complex<double> (0, 0)},
			    {complex<double> (0, 0), complex<double> (0, 0), complex<double> (0, 0), complex<double> (-1, 0)},
			    {complex<double> (-1, 0), complex<double> (0, 0), complex<double> (0, 0), complex<double> (0, 0)},
			    {complex<double> (0, 0), complex<double> (1, 0), complex<double> (0, 0), complex<double> (0, 0)} };

			    
const Matrix Id = Matrix(id);
const Matrix Gamma[4] = {Matrix(h0),Matrix(h1),Matrix(h2),Matrix(h3)};
const Matrix Gamma5 = I*Gamma[0]*Gamma[1]*Gamma[2]*Gamma[3];

const Matrix Gamma_mu5[4] = {Gamma[0]*Gamma5,
			  Gamma[1]*Gamma5,
			  Gamma[2]*Gamma5,
			  Gamma[3]*Gamma5};
			  
const Matrix Gamma_munu[4][4] = {{Gamma[0]*Gamma[0], Gamma[0]*Gamma[1], Gamma[0]*Gamma[2], Gamma[0]*Gamma[3] }, 
			  {Gamma[1]*Gamma[0], Gamma[1]*Gamma[1], Gamma[1]*Gamma[2], Gamma[1]*Gamma[3] },
			  {Gamma[2]*Gamma[0], Gamma[2]*Gamma[1], Gamma[2]*Gamma[2], Gamma[2]*Gamma[3] },
			  {Gamma[3]*Gamma[0], Gamma[3]*Gamma[1], Gamma[3]*Gamma[2], Gamma[3]*Gamma[3] },
			 };

const Matrix mGamma_munu[4][4] = {{(-1.)*Gamma[0]*Gamma[0], (-1.)*Gamma[0]*Gamma[1], (-1.)*Gamma[0]*Gamma[2], (-1.)*Gamma[0]*Gamma[3] }, 
			  {(-1.)*Gamma[1]*Gamma[0], (-1.)*Gamma[1]*Gamma[1], (-1.)*Gamma[1]*Gamma[2], (-1.)*Gamma[1]*Gamma[3] },
			  {(-1.)*Gamma[2]*Gamma[0], (-1.)*Gamma[2]*Gamma[1], (-1.)*Gamma[2]*Gamma[2], (-1.)*Gamma[2]*Gamma[3] },
			  {(-1.)*Gamma[3]*Gamma[0], (-1.)*Gamma[3]*Gamma[1], (-1.)*Gamma[3]*Gamma[2], (-1.)*Gamma[3]*Gamma[3] },
			 };
