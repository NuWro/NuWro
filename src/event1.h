#ifndef _event1_h_
#define _event1_h_

#include <iostream>
#include <vector>
#include "TObject.h"
#include "particle.h"
#include "params.h"

class flags
{
	public:
		/// primary vertex flags
		bool qel;        ///< (quasi) elastic       (qel == dyn/2==0)
		bool res;        ///< resontant             (res == dyn/2==1) 
		bool dis;        ///< deep inelastic        (dis == dyn/2==2)
		bool coh;        ///< coherent              (coh == dyn/2==3)
		bool mec;        ///< meson exhange current (mec == dyn/2==4)
		bool hyp;        ///< hyperon production    
    bool lep;        ///< neutrino-lepton                

		bool nc;         ///< neutral current       (nc == dyn%2)     
		bool cc;         ///< charged current       (cc == !nc)

		bool anty;       ///< true if antineutrino (anty==in[0].pdg<0)

		bool res_delta{false};  ///< true if RES pion comes from Delta
		bool need_resample_dir{false};
		bool need_resample_phi{false};
};

using namespace std;

class event:public TObject
{
	public:
		flags flag;               ///< flags for convenient filtering of the events in root scripts  
		params par;               ///< copy of all the input parameters 
		                          ///< NOTE: some parameters (e.g. nucleus_p, nucleus_n) can vary 		                          
		                          ///< from event to event (e.g. in case of the detector simulation)
		vector < particle > in;	  ///< vector of incoming particles
		vector < particle > temp; ///< vector of temporary particles (daughters of primary vertex in DIS)
		vector < particle > out;  ///< vector of outgoing particles (before fsi)
		vector < particle > post; ///< vector of particles leaving the nucleus
		vector < particle > all;  ///< vector of all particles (inclluding temporary fsi particles)

		double weight;	        ///< cross section of current event in cm^2 (set to total cross section on saving to file)
		double norm;            ///< norm of the initial neutrino (for weighted beams; not used) 
		vec r;                  ///< position of the event inside the detector
		double density;         ///< density of the detector matter in point of the interaction
		int dyn;		        ///< dynamics channel of the primary vertex. Possible values are:
		                        ///< 0,1 - qel  cc/nc - (quasi) elastic
		                        ///< 2,3 - res  cc/nc - resonant (via delta) with some background
		                        ///< 4,5 - dis  cc/nc - deep inelastric
		                        ///< 6,7 - coh  cc/nc - coherent 
		                        ///< 8,9 - mec  cc/nc - meson exhchange current 
		                        ///< 12 - lep  cc/nc - neutrino-lepton interaction
		                        ///< 20 -  eel  nc - elastic electron scattering 

		int nod[18];            ///< number of rescattering interactions of given type:
		                        ///< 0 - nucleon elastic,
		                        ///< 1 - nucleon ce,
		                        ///< 2 - nucleon spp,
		                        ///< 3 - nucleon dpp,
		                        ///< 4 - pion elastic,
		                        ///< 5 - pion ce,
		                        ///< 6 - pion spp,
		                        ///< 7 - pion dpp,
		                        ///< 8 - pion abs,
		                        ///< 9 - jailed,
		                        ///< 10 - escape
		                        ///< 11 - pion tpp
		                        ///< 12 - pion no interaction
		                        ///< 13 - nucleon no interaction
					///< 14 - hyperon no interaction
					///< 15 - hyperon elastic scatter
					///< 16 - hyperon Lambda -> Sigma conversion
					///< 17 - hyperon Sigma -> Lambda conversion
		int pr;     ///< number of protons  in the residual nucleus
		int nr;     ///< number of neutrons in the residual nucleus
		double r_distance; //< distance from nucleus center of absorption point (if happened)

		double res_jacobian; ///< store Jacobian calculated in RES for random kinematics
		double res_angrew;   ///< store xsec factor coming from angular distribution (for Delta)
		particle res_nu;     ///< store neutrino for reweighting
		vect res_q;          ///< store q for reweighting

		event ():weight(0),norm(1){}///< default constructor
		inline void check();        ///< stop program if event weight or momentum of any particle is NaN (not a number) 
		inline void clear_fsi();    ///< clear the fsi intermediate particles tracks
		inline particle nu ();	    ///< initial neutrino
		inline particle N0 ();	    ///< initial nucleon
		inline vect q ();		    ///< fourmomentum transfer
		inline double q0 ();	    ///< energy transfer
		inline double qv ();	    ///< momentum transfer
		inline double q2 ();	    ///< momentum transfer squared
		inline double s ();         ///< s - variable 
		inline double costheta ();  ///< cos theta lab
		inline double E ();         ///< initial neutrino energy
		inline int charge (int r);  ///< total charge: 0 - initial, 1 - before fsi, 2 - after fsi
		inline double W ();         ///< invariant mass (before fsi, all particles except the rescattered lepton)
		inline int n ();            ///< number of particles after primery vertex (before fsi)
		inline int f ();            ///< number of particles leaving nucleous
		inline int nof (int pdg);                     ///< number of particles after primery vertex 
		inline int nof (int pdg1,int pdg2);           ///< number of particles after primery vertex 
		inline int nof (int pdg1,int pdg2,int pdg3);  ///< number of particles after primery vertex 
		inline int fof (int pdg);                     ///< number of particles leaving nucleus
		inline int fof (int pdg1, int pdg2);          ///< number of particles leaving nucleus
		inline int fof (int pdg1, int pdg2, int pdg3);///< number of particles leaving nucleus
		inline double przod ();
		inline double tyl ();
		inline int number_of_nucleon_elastic ();      ///< number of nucleon elastic interactions during fsi
		inline int number_of_nucleon_ce ();           ///< number of nucleon ce interactions during fsi
		inline int number_of_nucleon_spp ();          ///< number of nucleon spp interactions during fsi
		inline int number_of_nucleon_dpp ();          ///< number of nucleon dpp interactions during fsi
		inline int number_of_pion_elastic ();         ///< number of pion elastic interactions during fsi
		inline int number_of_pion_ce ();              ///< number of pion ce interactions during fsi
		inline int number_of_pion_spp ();             ///< number of pion spp interactions during fsi
		inline int number_of_pion_dpp ();             ///< number of pion dpp interactions during fsi
		inline int number_of_pion_tpp ();             ///< number of pion tpp interactions during fsi
		inline int number_of_pion_abs ();             ///< number of pions absorbed during fsi
		inline int number_of_pion_no_interactions (); ///< number of pion steps with no interactions
		inline int number_of_nucleon_no_interactions (); ///< number of nucleon steps with no interactions
		inline double absorption_position ();         ///< positions where absorption occured
		inline int number_of_jailed ();               ///< number of nucleons jailed in nucleous during fsi
		inline int number_of_escape ();               ///< number of particles that escaped from nucleus during fsi
		inline int number_of_interactions ();         ///< total number of interactions during fsi
		inline int number_of_particles (int pdg, bool fsi); ///< number of particles before/after fsi
		inline double nuc_kin_en();                         ///< total kinetic energy of nucleons that left the nucleus
		inline int num_part_thr (int pdg, bool fsi, double threshold);  ///< number of particles with momentum above threshold before/after fsi
		inline int num_part_thr_withincosine (int pdg, bool fsi, double threshold, double kosinus);/// as above plus a condition for cosine
		inline int num_part_two_thr_withincosine (int pdg, bool fsi, double threshold_min, double theshold_max, double kosinus);/// as above but minimum and maximum momentum thresholds
		inline double proton_cosine(bool fsi, double thr);              ///< cosine of the angle between two protons with momenta above thr
		inline double proton_transp_mom();
		inline double proton_transp_mom2();
		inline int proton_transp();
		inline int proton_pair_number1 (bool fsi, double thr);
		inline int proton_pair_number2 (bool fsi, double thr);
		inline double part_max_mom (int pdg, bool fsi);         ///< maximal momentum of particle pdg
		inline double part_sec_mom (int pdg, bool fsi);         ///< second largest momentum of particle pdg 
		inline double vert_act (double pion_threshold, bool fsi, double proton_threshold);
		inline double Erec (double Bin);                        ///< reconstructed neutrino energy 
		inline double Q2rec (double Bin); ///reconstructed Q2
		inline double proton_recoil ();
		inline double neutron_recoil ();
		inline double photon_recoil ();
		inline double meson_recoil_without_masses ();
		inline double meson_recoil_with_masses ();
		inline double lepton_recoil ();
		inline double total_recoil_with_masses (double K, double N);
		inline double total_recoil_without_masses (double K, double N);
		inline double neutral_kaon_recoil ();
		inline vec proton_max_mom();
		inline vect particle_max_mom(int pdg, bool fsi);
		inline vect particle_max_mom_withincosine (int pdg, bool fsi, double kosinus);
		inline vect particle_max_mom_withincosine_withinmomentum (int pdg, bool fsi, double kosinus, double threshold_min, double threshold_max);
		inline double total_hadr_post();
		ClassDef (event, 1);
};

/// add particle to a list and set its id to position in the list.
inline void registration (vector<particle>& list, particle &p);

/// I M P L E M E N E T A T I O N

/// clear the fsi intermediate particles tracks
void event::clear_fsi()
{
	post=vector<particle>();
	all=vector<particle>();
	//for(int i=0;i<12;i++) nod[i]=0; //it is already done on the begining of the cascade
}


/// initial neutrino
particle event::nu ()
{
	return in[0];
}


/// initial nucleon
particle event::N0 ()
{
	return in[1];
}


/// fourmomentum transfer
vect event::q ()
{
	vect q = in[0] - out[0];
	return q;
}


/// energy transfer
double event::q0 ()
{
	return in[0].t - out[0].t;
}


/// momentum transfer
double event::qv ()
{
	vect q = in[0] - out[0];
	return vec(q).norm();
}


/// fourmomentum transfer squared
double event::q2 ()
{
	vect q = in[0] - out[0];
	return q * q;
}


/// s variable
double event::s ()
{
	vect q = in[0] + in[1];
	return q * q;
}


/// cos theta lab
double event::costheta ()
{
	return cos(vec(in[0]),vec(out[0]));
}


/// neutrino energy 
double event::E ()
{
	return in[0].t;
}


/// number of particles after primary vertex
int event::n ()
{
	return out.size ();
}


/// number of particles leaving nucleus
int event::f ()
{
	return post.size ();
}


/// number of particles with given pdg after primary vertex 
int event::nof (int pdg)
{
	int c = 0;
	for (int i = 0; i < out.size (); i++)
		if(out[i].pdg == pdg)
			c++;
	return c;
}


/// number of particles with given pdg after primary vertex 
int event::nof (int pdg1,int pdg2)
{
	int c = 0;
	for (int i = 0; i < out.size (); i++)
		if(out[i].pdg == pdg1 || out[i].pdg == pdg2)
			c++;
	return c;
}


/// number of particles with given pdg after primary vertex 
int event::nof (int pdg1,int pdg2,int pdg3)
{
	int c = 0;
	for (int i = 0; i < out.size (); i++)
		if (out[i].pdg == pdg1 || out[i].pdg == pdg2 || out[i].pdg == pdg3)
			c++;
	return c;
}


/// number of particles with given pdg leaving nucleus
int event::fof (int pdg)
{
	int c = 0;
	for (int i = 0; i < post.size (); i++)
		if(post[i].pdg == pdg)
			c++;
	return c;
}


/// number of particles with given pdg leaving nucleus
int event::fof (int pdg1, int pdg2)
{
	int c = 0;
	for (int i = 0; i < post.size (); i++)
		if(post[i].pdg==pdg1 || post[i].pdg==pdg2)
			c++;
	return c;
}


/// number of particles with given pdg leaving nucleus
int event::fof (int pdg1, int pdg2, int pdg3)
{
	int c = 0;
	for (int i = 0; i < post.size (); i++)
		if(post[i].pdg==pdg1 || post[i].pdg==pdg2 || post[i].pdg==pdg3)
			c++;
	return c;
}

/// total charge: 0 - initial, 1 - before fsi, 2 - after fsi
int event::charge (int r)
{
	vector<particle>  *set=&post;
	int c = 0;
	switch(r)
	{
		case 0: return par.nucleus_p+in[0].charge();break;
		case 1: set=&out;c=par.nucleus_p-(dyn<6)*in[1].charge();break;
		case 2: set=&post;c=pr;break;
	}
	for (int i = 0; i < set->size (); i++)
		c+=(*set)[i].charge();
	return c;
}

/// invaraint mass 
double event::W ()
{
	vect h = out[1];
	for (int a = 2; a < out.size (); a++)
		h = h + out[a];
	return sqrt (h * h);
}

/// 
double event::przod ()
{
	int licz = 0;
	vect h = out[1];
	for (int a = 2; a < out.size (); a++)
		h = h + out[a];

	vect tran = in[0] - out[0];
	vec ptr = vec (tran.x, tran.y, tran.z);

	for (int b = 1; b < out.size (); b++)
	{
		if ((out[b].pdg != 111) && (out[b].pdg != 22) && (out[b].pdg != 130)
			&& (out[b].pdg != 2112) && (out[b].pdg != -2112))

		{
			vect hh = out[b];
			hh.boost (-h.v());

			vec phh = vec (hh.x, hh.y, hh.z);

			if (ptr * phh > 0)
				licz = licz + 1;
		}
	}
	return licz;
}

/// 
double event::tyl ()
{
	int licz = 0;
	vect h = out[1];
	for (int a = 2; a < out.size (); a++)
		h = h + out[a];

	vect tran = in[0] - out[0];
	vec ptr = vec (tran.x, tran.y, tran.z);

	for (int b = 1; b < out.size (); b++)
	{
		if ((out[b].pdg != 111) && (out[b].pdg != 22) && (out[b].pdg != 130)
			&& (out[b].pdg != 2112) && (out[b].pdg != -2112))

		{
			vect hh = out[b];
			hh.boost (-(h.v ()));

			vec phh = vec (hh.x, hh.y, hh.z);

			if (ptr * phh < 0)
				licz = licz + 1;
		}
	}
	return licz;
}


int event::number_of_nucleon_elastic ()
{
	return nod[0];
}


int event::number_of_nucleon_ce ()
{
	return nod[1];
}


int event::number_of_nucleon_spp ()
{
	return nod[2];
}


int event::number_of_nucleon_dpp ()
{
	return nod[3];
}


int event::number_of_pion_elastic ()
{
	return nod[4];
}


int event::number_of_pion_ce ()
{
	return nod[5];
}


int event::number_of_pion_spp ()
{
	return nod[6];
}


int event::number_of_pion_dpp ()
{
	return nod[7];
}


int event::number_of_pion_tpp ()
{
	return nod[11];
}


int event::number_of_pion_abs ()
{
	return nod[8];
}


int event::number_of_jailed ()
{
	return nod[9];
}


int event::number_of_escape ()
{
	return nod[10];
}

int event::number_of_pion_no_interactions ()
{
	return nod[12];
}

int event::number_of_nucleon_no_interactions ()
{
	return nod[13];
}

double event:: absorption_position ()
{
  return r_distance;
}


int event::number_of_interactions ()
{
	int noi = 0;
	for (int i = 0; i<11; i++)
	{
		noi = noi + nod[i];
	}

	return noi;
}


/// number of particles of given pdg code ( 0 - before FSI, 1 - after FSI )
int event::number_of_particles (int pdg, bool fsi)
{
	int number = 0;

	if(fsi)
	{
		for (int k = 0; k<post.size(); k++)
		{
			if (post[k].pdg == pdg) number++;
		}
	}
	else if (!fsi)
	{
		for (int k = 0; k<out.size(); k++)
		{
			if (out[k].pdg == pdg) number++;
		}
	}

	return number;

}

/// total kinetic energy of nucleons that left the nucleus
double event::nuc_kin_en()
{
	double sum=0;
	for (int k = 0; k<post.size(); k++)
	{
		if ( post[k].pdg == 2212 || post[k].pdg == 2112 )
			sum+=post[k].Ek();
	}
	return sum;
}

double event::proton_recoil()//contain very small antineutron contribution
{
  double sum=0;
	for (int k = 0; k<post.size(); k++)
	{
		if ( post[k].pdg == 2212 )
			sum+=post[k].Ek();
		if ( post[k].pdg == -2212 )
			sum+=post[k].Ek()+2*post[k].mass();
		if ( post[k].pdg == -2112 )
			sum+=post[k].Ek()+2*post[k].mass();
	}
	return sum;
}

double event::neutron_recoil()
{
  double sum=0;
	for (int k = 0; k<post.size(); k++)
	{
		if ( post[k].pdg == 2112 )
			sum+=post[k].Ek();
		
	}
	return sum;
}

double event::photon_recoil()
{
  double sum=0;
	for (int k = 0; k<post.size(); k++)
	{
		if ( post[k].pdg == 22 )
			sum+=post[k].t;
	}
	return sum;
}

double event::lepton_recoil()
{
  double sum=0;
  if (post.size()>0)
  {
	for (int k = 1; k<post.size(); k++)
	{
		if ( post[k].pdg == 13 || post[k].pdg == -13 || post[k].pdg == 11 || post[k].pdg == -11 )
			sum+=post[k].t;
	}
  }
	return sum;
}

double event::meson_recoil_with_masses()
{
  double sum=0;
	for (int k = 0; k<post.size(); k++)
	{
		if ( post[k].pdg == 211 || post[k].pdg == -211 ) 
			sum+=post[k].t;
		if ( post[k].pdg == 111 || post[k].pdg == 321 || post[k].pdg == -321 )
			sum+=post[k].t;
	}
	return sum;
}

double event::meson_recoil_without_masses()
{
  double sum=0;
	for (int k = 0; k<post.size(); k++)
	{
		if ( post[k].pdg == 211 || post[k].pdg == -211 ) 
			sum+=post[k].Ek();
		if ( post[k].pdg == 111 || post[k].pdg == 321 || post[k].pdg == -321 )
			sum+=post[k].t;
	}
	return sum;
}

double event::neutral_kaon_recoil()
{
  double sum=0;
	for (int k = 0; k<post.size(); k++)
	{
		if ( post[k].pdg == 311 || post[k].pdg == -311  || post[k].pdg ==130 || post[k].pdg == 310 )
			sum+=post[k].t;
	}
	return sum;
}


double event::total_recoil_with_masses (double K0_fraction, double neutron_fraction)
{
  return meson_recoil_with_masses () + lepton_recoil() + photon_recoil() + proton_recoil() + 
  neutron_fraction*neutron_recoil() + K0_fraction*neutral_kaon_recoil() ;
}

double event::total_recoil_without_masses (double K0_fraction, double neutron_fraction)
{
  return meson_recoil_without_masses () + lepton_recoil() + photon_recoil() + proton_recoil() + 
  neutron_fraction*neutron_recoil() + K0_fraction*neutral_kaon_recoil() ;
}

double event::total_hadr_post()
{
  double wynik=0;
	for (int a = 1; a < post.size (); a++)
		wynik+=post[a].t;
	return wynik;
}


/// number of particles of given pdg with momentum above threshold before/after fsi 
int event::num_part_thr (int pdg, bool fsi, double threshold)
{
	int number = 0;

	if(fsi)
	{
		for (int k = 0; k<post.size(); k++)
		{
			if ( (post[k].pdg == pdg) && (post[k].momentum() > threshold) )
				number++;
		}
	}
	else if (!fsi)
	{
		for (int k = 0; k<out.size(); k++)
		{
			if ( (out[k].pdg == pdg) && (out[k].momentum() > threshold) )
				number++;
		}
	}

	return number;

}

/// number of particles of given pdg with momentum above threshold before/after fsi and with angle within cosine region
int event::num_part_thr_withincosine (int pdg, bool fsi, double threshold, double kosinus)
{
	int number = 0;

	if(fsi)
	{
		for (int k = 0; k<post.size(); k++)
		{
			if ( post[k].pdg == pdg && post[k].momentum() > threshold && post[k].z/post[k].momentum() > kosinus )
				number++;
		}
	}
	else if (!fsi)
	{
		for (int k = 0; k<out.size(); k++)
		{
			if ( out[k].pdg == pdg && out[k].momentum() > threshold && out[k].z/out[k].momentum() > kosinus )
				number++;
		}
	}

	return number;
}

/// number of particles of given pdg with momentum above threshold before/after fsi and with angle within cosine region
int event::num_part_two_thr_withincosine (int pdg, bool fsi, double threshold_min, double threshold_max, double kosinus)
{
	int number = 0;

	if(fsi)
	{
		for (int k = 0; k<post.size(); k++)
		{
			if ( post[k].pdg == pdg && post[k].momentum() > threshold_min && post[k].momentum() < threshold_max && post[k].z/post[k].momentum() > kosinus )
				number++;
		}
	}
	else if (!fsi)
	{
		for (int k = 0; k<out.size(); k++)
		{
			if ( out[k].pdg == pdg && out[k].momentum() > threshold_min && out[k].momentum() < threshold_max && out[k].z/out[k].momentum() > kosinus )
				number++;
		}
	}

	return number;
}



/// cosine if the angle between the two outgoing protons
double event::proton_cosine(bool fsi, double thr)
{
	int numer[2];
	int ile=0;
	if (fsi)
	{
		for (int k = 0; k<post.size(); k++)
		{
			if ( post[k].pdg == 2212 && post[k].momentum() > thr  )
			{
				numer[ile]=k;
				ile++;
			}
		}
		return cos(post[numer[0]],post[numer[1]]);
		// return ( post[numer[0]].x*post[numer[1]].x + post[numer[0]].y*post[numer[1]].y + post[numer[0]].z*post[numer[1]].z )/post[numer[0]].momentum()/post[numer[1]].momentum();

	}
	else
	{
		for (int k = 0; k<out.size(); k++)
		{
			if ( out[k].pdg == 2212 && out[k].momentum() > thr  )
			{
				numer[ile]=k;
				ile++;
			}
		}
		return cos(out[numer[0]],out[numer[1]]);
		//return ( out[numer[0]].x*out[numer[1]].x + out[numer[0]].y*out[numer[1]].y + out[numer[0]].z*out[numer[1]].z )/out[numer[0]].momentum()/out[numer[1]].momentum();

	}

}


///calculates momentum of proton which did not suffer from fsi (the first one, sometimes there are two,see below)
double event::proton_transp_mom()
{
	for (int k = 0; k<out.size(); k++)
		if (out[k].pdg ==2212)
	{
		for (int l = 0; l<post.size(); l++)
		{
			if (post[l].pdg ==2212)
			{
//				double kos = ( out[k].x*post[l].x + out[k].y*post[l].y + out[k].z*post[l].z )/out[k].momentum()/post[l].momentum();
				double kos = cos(out[k],post[l]);
				if (kos>0.999)
					return out[k].momentum();
			}
		}
	}
	return 0;
}


///calculates momentum of proton which did not suffer from fsi (the second one if it happens -- very rarely, but stil...)
double event::proton_transp_mom2()
{
	int counter=0;
	for (int k = 0; k<out.size(); k++)
		if (out[k].pdg ==2212)
	{
		for (int l = 0; l<post.size(); l++)
		{
			if (post[l].pdg ==2212)
			{
				double kos = ( out[k].x*post[l].x + out[k].y*post[l].y + out[k].z*post[l].z )/out[k].momentum()/post[l].momentum();
				if (kos>0.999 && counter==1)
					return out[k].momentum();

				if (kos>0.999 && counter==0)
					counter++;
			}
		}
	}
	return 0;
}


///calculates number of protons which did not suffer from fsi
int event::proton_transp()
{
	int ile =0;
	for (int k = 0; k<out.size(); k++)
		if (out[k].pdg ==2212)
		{
			for (int l = 0; l<post.size(); l++)
			{
				if (post[l].pdg ==2212)
				{
//					double kos = ( out[k].x*post[l].x + out[k].y*post[l].y + out[k].z*post[l].z )/out[k].momentum()/post[l].momentum();
					double kos = cos(out[k],post[l]);
					if (kos>0.999)
						ile++;
				}
			}
		}
	return ile;
}


int event::proton_pair_number1 (bool fsi, double thr)
{
	int numer[2];
	int ile=0;
	if (fsi)
	{
		for (int k = 0; k<post.size(); k++)
		{
			if ( post[k].pdg == 2212 && post[k].momentum() > thr  )
			{
				numer[ile]=k;
				ile++;
			}
		}
		return numer[0];
	}
	else
	{
		for (int k = 0; k<out.size(); k++)
		{
			if ( out[k].pdg == 2212 && out[k].momentum() > thr  )
			{
				numer[ile]=k;
				ile++;
			}
		}
		return numer[0];
	}
}


int event::proton_pair_number2 (bool fsi, double thr)
{
	int numer[2];
	int ile=0;
	if (fsi)
	{
		for (int k = 0; k<post.size(); k++)
		{
			if ( post[k].pdg == 2212 && post[k].momentum() > thr  )
			{
				numer[ile]=k;
				ile++;
			}
		}
		return numer[1];
	}
	else
	{
		for (int k = 0; k<out.size(); k++)
		{
			if ( out[k].pdg == 2212 && out[k].momentum() > thr  )
			{
				numer[ile]=k;
				ile++;
			}
		}
		return numer[1];
	}
}


/// maximal momentum of particle with given pdg before/after FSI (when fsi=0/1 resp.)
double event::part_max_mom (int pdg, bool fsi)
{
	double mom = 0.0;

	if(fsi)
	{

		for (int k = 0; k<post.size(); k++)
		{
			if ( (post[k].pdg == pdg) && (post[k].momentum() > mom) )
				mom=post[k].momentum();
		}
	}
	else if (!fsi)
	{
		for (int k = 0; k<out.size(); k++)
		{
			if ( (out[k].pdg == pdg) && (out[k].momentum() > mom) )
				mom=out[k].momentum();
			//cout<<mom<<endl;
		}
	}

	return mom;

}

vec event::proton_max_mom()
{
  vec proton_mom;
  double length=0;

	for (int k = 0; k<post.size(); k++)
	{
	if ( (post[k].pdg == 2212) && (post[k].momentum() > length) )
	{proton_mom=post[k].p();
	  length = post[k].momentum();
	//cout<<"sabaka  "<<proton_mom<<"  "<<length<<endl;
	}
	}
	return proton_mom;
}

vect event::particle_max_mom(int pdg, bool fsi)
{
  vect particle_mom;
  double length=0;
    if (fsi==true)
    {
	for (int k = 0; k<post.size(); k++)
	{
	if ( (post[k].pdg == pdg) && (post[k].momentum() > length) )
	{particle_mom=post[k];
	  length = post[k].momentum();
	//cout<<"sabaka  "<<proton_mom<<"  "<<length<<endl;
	}
	}
    }
	if (fsi==false)
    {
	for (int k = 0; k<out.size(); k++)
	{
	if ( (out[k].pdg == pdg) && (out[k].momentum() > length) )
	{particle_mom=out[k];
	  length = out[k].momentum();
	//cout<<"sabaka  "<<proton_mom<<"  "<<length<<endl;
	}
	}
    }
	
	return particle_mom;
}


vect event::particle_max_mom_withincosine (int pdg, bool fsi, double kosinus)
{
  vect particle_mom;
  double length=0;
    if (fsi==true)
    {
	for (int k = 0; k<post.size(); k++)
	{
	if ( post[k].pdg == pdg && post[k].momentum() > length && post[k].z/post[k].momentum() > kosinus)
	{particle_mom=post[k];
	  length = post[k].momentum();
	//cout<<"sabaka  "<<proton_mom<<"  "<<length<<endl;
	}
	}
    }
	if (fsi==false)
    {
	for (int k = 0; k<out.size(); k++)
	{
	if ( out[k].pdg == pdg && out[k].momentum() > length && out[k].z/out[k].momentum() > kosinus )
	{particle_mom=out[k];
	  length = out[k].momentum();
	//cout<<"sabaka  "<<proton_mom<<"  "<<length<<endl;
	}
	}
    }
	
	return particle_mom;
}

vect event::particle_max_mom_withincosine_withinmomentum (int pdg, bool fsi, double kosinus, double threshold_min, double threshold_max)
{
  vect particle_mom;
  double length=0;
    if (fsi==true)
    {
	for (int k = 0; k<post.size(); k++)
	{
	if ( post[k].pdg == pdg && post[k].momentum() > threshold_min && post[k].momentum() < threshold_max && post[k].momentum() > length && post[k].z/post[k].momentum() > kosinus)
	{particle_mom=post[k];
	  length = post[k].momentum();
	//cout<<"sabaka  "<<proton_mom<<"  "<<length<<endl;
	}
	}
    }
	if (fsi==false)
    {
	for (int k = 0; k<out.size(); k++)
	{
	if ( out[k].pdg == pdg && out[k].momentum() > threshold_min && out[k].momentum() < threshold_max && out[k].momentum() > length && out[k].z/out[k].momentum() > kosinus )
	{particle_mom=out[k];
	  length = out[k].momentum();
	//cout<<"sabaka  "<<proton_mom<<"  "<<length<<endl;
	}
	}
    }
	
	return particle_mom;
}


/// second largest momentum of particle with given pdg before/after FSI (when fsi=0/1 resp.)
double event::part_sec_mom (int pdg, bool fsi) 
{
	double mom [2]= {0.0, 0.0};
	double memory;

	if(fsi)
	{
		for (int k = 0; k<post.size(); k++)
		{
			if  (post[k].pdg == pdg)
			{
				double ped =post[k].momentum();
				if (ped >mom[1])
					{mom[1]=ped;}
					if (ped>mom[0])
				{
					memory = mom[0];
					mom[0]=ped;
					mom[1]=memory;
				}
			}
		}
	}
	else if (!fsi)
	{
		for (int k = 0; k<out.size(); k++)
		{
			if  (out[k].pdg == pdg)
			{
				double ped =out[k].momentum();
				//cout<<ped<<endl;
				if (ped >mom[1])
					{mom[1]=ped;}
					if (ped>mom[0])
				{
					memory = mom[0];
					mom[0]=ped;
					mom[1]=memory;
				}
			}
		}
	}

	return mom[1];

}


double event::vert_act (double pion_threshold, bool fsi, double proton_threshold)
{
	double veract=0.0;

	if(fsi)
	{
		for (int k = 0; k<post.size(); k++)
		{
			if ((post[k].pdg == 211) && (post[k].momentum() < pion_threshold) )
				veract+=post[k].t-post[k].mass();

			if ((post[k].pdg == -211) && (post[k].momentum() < pion_threshold) )
				veract+=post[k].t-post[k].mass();

			if ((post[k].pdg == 2212) && (post[k].momentum() < proton_threshold) )
				veract+=post[k].t-post[k].mass();

								 //positive kaon; assume never reconstructed
			if (post[k].pdg == 321)
				veract+=post[k].t-post[k].mass();
		}
		return veract;
	}
	else
	{
		for (int k = 0; k<out.size(); k++)
		{
			if ((out[k].pdg == 211) && (out[k].momentum() < pion_threshold) )
				veract+=out[k].t-out[k].mass();

			if ((out[k].pdg == -211) && (out[k].momentum() < pion_threshold) )
				veract+=out[k].t-out[k].mass();

			if ((out[k].pdg == 2212) && (out[k].momentum() < proton_threshold) )
				veract+=out[k].t-out[k].mass();

								 //positive kaon; assume never reconstructed
			if (out[k].pdg == 321)
				veract+=out[k].t-out[k].mass();
		}
		return veract;
	}

}

/// Reconstructed neutrino energy 
double event::Erec (double Bin)
{
	double massprim = in[1].mass() - Bin;
	return ( out[0].t*massprim + 0.5* (out[1].mass()*out[1].mass() - out[0].mass()*out[0].mass() - massprim*massprim) )/
		( massprim - out[0].t + out[0].z);
}

double event::Q2rec (double Bin)
{
	return -out[0].mass()*out[0].mass() + 2.0*Erec(Bin)*(out[0].t -out[0].z);
}


/// stop program if event weight or momentum of any particle is NaN (not a number) 
void event::check()
{
	for(int i=0;i<out.size();i++)
	{
		particle p=out[i];
		if(!(p.x==p.x && p.y==p.y && p.z==p.z and p.t==p.t) )
		{
			cerr<<p<<endl;
			exit(20);
		}
	}
	if(!(weight==weight))
	{
		cerr<<"dyn="<<dyn<<endl;
		exit(21);
	}

}

/// add particle to a list and set its id to position in the list.
void registration (vector<particle>& list, particle &p)
{
	p.id=list.size();
	list.push_back(p);
}


#endif
