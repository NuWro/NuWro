#ifndef _LeptonMass_h_
#define _LeptonMass_h_
#include "vec.h"
#include "vect.h"

extern double lepton_mass(int lepton_in, bool cur);

void   kinfinder(vec b, vec& d, double kos);

double kin1part(double hama, int nukleon, int meson, vect& finnuk, vect& finpion, vec kierunek);

void   kin2part(double hama, int nukleon, int meson, vect& finnuk, vect& finpion);

void   kin3part(vect neutr, vect finlep, double hama, int nukleon2, int meson, vect& finnuk, vect& finpion);

double kin4part(vect neutr, vect finlep, double hama, int nukleon2, int meson, vect& finnuk, vect& finpion, int ANLang);

void   rotation(vect& cztero, vec trzy);

double binen(vec mom, int p, int n);

double binen2(double pp, int p, int n);

double deuter_binen(vec mom);

double deuter_binen2(double pp);

#endif
