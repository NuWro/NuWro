#ifndef _target_mixer_h_
#define _target_mixer_h_
#include "params.h"

double get_Eb(int p, int n);
double get_kf(int p, int n);
double get_model(int p, int n);

class target_mixer
{
protected:
   int N;
   int p[10];
   int n[10];
   
   int model[10];
   int Eb[10];
   int kf[10];
   double w[10];
   double sum;
//   params &par;
public:   
   target_mixer(params& par):N(0),sum(0)
   { stringstream s(par.target_content);
     while (s && N<10)
     { int p1,n1;
       double ratio;
       char perc;
       string line;
       getline(s,line);
       stringstream s1(line);
       s1>> p1>>n1 >> ratio>> perc;
       if(!s1) return;
       if(perc!='x') 
         {cerr<< "target_content= "<<p1<<' '<<n1<<' '<<ratio<<" 'x' sign expected"<<endl;
         exit(32);
         }
        Eb[N]=kf[N]=model[N]=0;
        s1>> Eb[N]>>kf[N]>>model[N];
        p[N]=p1;
        n[N]=n1;
       sum+=w[N]=ratio*(p1+n1);
       N++;
     } 
   }
      
   /// choose material for next event by adjusting params
   void prepare(params& par)
   { if(N==0) return;
	 double x=frandom()*sum,x1=0;
     int i=0;
     while((x1+=w[i])<x) i++;
     //cout<<" prepare "<<p[i] <<" "<<n[i]<<endl;
     par.nucleus_p=p[i];
     par.nucleus_n=n[i];
//     par.nucleus_E_b= Eb[i] ? Eb[i] : get_Eb(p[i],n[i]);
//     par.nucleus_kf= kf[i] ? kf[i] : get_kf(p[i],n[i]);
//     par.nucleus_target=model[i] ? model[i] : get_model(p[i],n[i]);
     par.nucleus_E_b= Eb[i]? Eb[i]:get_Eb(p[i],n[i]);
     par.nucleus_kf=  kf[i]?kf[i]:get_kf(p[i],n[i]);
    //  par.nucleus_target=model[i]?model[i]:get_model(p[i],n[i]);
     
//     par.list();
   }
   
};
#endif 
