#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>
using namespace std;
#include "singlepion.h"

double
SPPF (int j, int k, int l, int t, double W)
{
  if (W < SPP_MIN)
    return SPP[j][k][l][t][0];

  if (W >= SPP_MAX)
    {
      cout << "warning! W out of allowed range!" << endl;
      return SPP[j][k][l][t][NSPPbins-1];
    }

  if (W > SPP_MIN && W < SPP_MAX)
    {
      double Wdiff = W - SPP_MIN;
      int bin = int (Wdiff / NSPPSize);
      double surplus = Wdiff - bin * NSPPSize;
      return (1 - surplus / NSPPSize) * SPP[j][k][l][t][bin] +
	surplus / NSPPSize * SPP[j][k][l][t][bin + 1];
    }
   return 0;  
}


double
alfadis (int j, int k, int l, int t, double W, double blending_start, double blending_end)
{
  double alfa = 0;
  double & W_min = blending_start;
  double & W_max = blending_end;
//cout<<j<<"  "<<k<<"  "<<l<<endl;
  if (k == 1 || ((j == 0 && l == 0) || (j == 1 && l == 1)))
    {				//cout<<"a"; 
      alfa = 0.0;
    }
  else
    {
      if ((k == 0)
	  && ((j == 0 && l == 1 && t == 1) || (j == 1 && l == 0 && t == 1)))
	{			//cout<<"b"; 
	  alfa = 0.3;
	}
      else
	{
	  if ((k == 0)
	      && ((j == 0 && l == 1 && t == 0)
		  || (j == 1 && l == 0 && t == 2)))
	    {			//cout<<"c"; 
	      alfa = 0.2;
	    }
	}
    }

  if ((W > 1080) && W < W_min)
    {				//cout<<alfa<<endl; 
      return alfa * (W - 1080) / (W_min - 1080);
    }
  else
    {
      if (W >= W_min && W < W_max)
	{			//cout<<alfa<<endl; 
	  return alfa + (1 - alfa) * (W - W_min) / (W_max - W_min);
	}
      else
	{
	  if (W >= W_max)
	    {			//cout<<alfa<<endl; 
	      return 1;
	    }
	}
    }
    
    return 0;
}

double betadis (int j, int k, int l, int t, double W, double bkgr, double blending_start, double blending_end)
{
  double alfa = 0;
  double & W_min = blending_start;
  double & W_max = blending_end;
//cout<<j<<"  "<<k<<"  "<<l<<endl;
  if (k == 1 || ((j == 0 && l == 0) || (j == 1 && l == 1)))
    {				//cout<<"a"; 
      alfa = 0.0;
      W_min = 1300 - 75*bkgr;
    }
  else
    {
      if ((k == 0)
	  && ((j == 0 && l == 1 && t == 1) || (j == 1 && l == 0 && t == 1)))
	{			//cout<<"b"; 
	  alfa = 0.3 + 0.15*bkgr;
	}
      else
	{
	  if ((k == 0)
	      && ((j == 0 && l == 1 && t == 0)
		  || (j == 1 && l == 0 && t == 2)))
	    {			//cout<<"c"; 
	      alfa = 0.2 + 0.15*bkgr;
	    }
	}
    }

  if ((W > 1080) && W < W_min)
    {				//cout<<alfa<<endl; 
      return alfa * (W - 1080) / (W_min - 1080);
    }
  else
    {
      if (W >= W_min && W < W_max)
	{			//cout<<alfa<<endl; 
	  return alfa + (1 - alfa) * (W - W_min) / (W_max - W_min);
	}
      else
	{
	  if (W >= W_max)
	    {			//cout<<alfa<<endl; 
	      return 1;
	    }
	}
    }
    
    return 0;
}


double
alfadelta (int j, int k, int l, int t, double W, double blending_start, double blending_end)
{
  return 1 - alfadis (j, k, l, t, W, blending_start, blending_end);
}
