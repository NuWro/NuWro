#ifndef _EnergyProfile_h_
#define _EnergyProfile_h_
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
#include "generatormt.h"
#include "jednostki.h"
#include "TFile.h"
#include "TH1.h"
#include "dirs.h"

using namespace std;

////////////////////////////////////////////////////////////////////////////////
/// EnergyProfile shoots a double precision number whose probability density is
/// given by the parameters to the constructor
////////////////////////////////////////////////////////////////////////////////

class EnergyProfile {
 private:
  int n;
  double Emin;
  double Emax;
  std::vector<double> prob, prob2;

  TH1D* spectrum;
  std::vector<double> Elow, Ehigh;
  bool roothist;
  
 public:
  ~EnergyProfile() {
  }

  /// Beam of Zero Energy
 EnergyProfile() : n(0), Emin(0), Emax(0), prob({}), prob2({}), roothist(false)  {}

  /// Monoenergetic beam
 EnergyProfile(double E) : n(0), Emin(E), Emax(E), prob({}), prob2({}), roothist(false) {}

  /// uniform probability density betweeen E1 and E2
  EnergyProfile(double E1, double E2)
    : n(1), Emin(E1), Emax(E2), prob({}), prob2({}), roothist(false) {}

  /// energy profile given by a histogram encoded in a string
  /// the string contains space separated values:
  /// Emin Emax b1 b2 ... bn
  /// where b1 b2 ... bn are the heightsof the bars of the histogram
  /// the widths are assumed equal to (Emax-Emin)/n
 EnergyProfile(string s) : n(0), Emin(0), Emax(0), prob({}), prob2({}), roothist(false) {
    read(s);
  }

  // Energy Profile From File and histogram name
 EnergyProfile(string f, string h) : n(0), Emin(0), Emax(0), prob({}), prob2({}), roothist(true) {
    readHistRoot(f, h, false);
  }

  // Set the Energy Profile from a ROOT hist
  void readHistRoot(string f, string h, bool inMEV);
  
  /// the auxiliary helper function for parsing the string paramater to
  /// the EnergyProfile constructor
  void read(string s);

  /// a random number conforming to the probability profile
  /// specified by construcor parameters
  double shoot(bool dis);

  /// diagnostic function for verifying if the beam
  /// has been correctly constructed
  void print();

  double minE();

  double GetN() const { return n; }
  double GetEmin() const { return Emin; }
  double GetEmax() const { return Emax; }
  double GetBinProb(int bin, bool dis) const {
    if (bin >= n) {
      std::cerr << "[ERROR]: Requested probability in bin: " << bin
                << " but this EnergyProfile only contains " << n << " bins."
                << std::endl;
      throw;
    }
    return (dis ? (bin ? (prob2[bin] - prob2[bin - 1]) : prob2[bin])
                : (bin ? (prob[bin] - prob[bin - 1]) : prob[bin]));
  }

  /// how many times weighting by E enhances this beam
  double disratio() {
    if (n == 0)
      return 1;
    else
      return prob2[n - 1] / prob[n - 1];
  }

  TH1D GetHist(bool inMEV=false) const;

};  // Energy Profile class

////////////////////////////////////////////////////////////////////////////////
///              I M P L E M E N T A I O N
////////////////////////////////////////////////////////////////////////////////

/// parse the string parameter to the EnergyProfile constructor
inline void EnergyProfile::read(string s) {
  prob.clear();
  prob2.clear();
  Emin = Emax = n = 0;

  stringstream in(s);
  in >> Emin >> Emax;
  Emin *= MeV;
  Emax *= MeV;
  if (!in || Emax==Emin) return;
  
  std::vector <double> bufor;
    
  while (in >> bufor.emplace_back()) 
    n++;

  if (n == 0) 
  {  
    bufor.resize(1);
    n = 1;
    bufor[0] = 1;
  }

  prob.resize(n);
  prob2.resize(n);
  Elow.resize(n);
  Ehigh.resize(n);
  
  double prev1 = 0, prev2 = 0;
  for (int i = 0; i < n; i++) {

    Elow[i] = Emin + (i) * (Emax - Emin)/n;
    Ehigh[i] = Emin + (i+1) * (Emax - Emin)/n;
    double E = (Elow[i] + Ehigh[i]) / 2.0;
  
    prob[i] = prev1 += bufor[i];
    prob2[i] = prev2 += bufor[i] * E;
  }
}

// Read in a flux histogram for the spectrum
inline void EnergyProfile::readHistRoot(string f, string h, bool inMEV){

  // Initial Setup
  prob.clear();
  prob2.clear();
  Emin = Emax = n = 0;

  // Read in histogram
  cout << "Opening TFile = " << f << endl;
  TFile* infile = new TFile(f.c_str(),"READ");
  if (infile->IsZombie()){
    cout << "Can't find flux file in: " << f << endl;
    cout << "Checking: " << get_data_dir() + "/beam/" + f << endl;

    infile = new TFile( (get_data_dir() + "/beam/" + f).c_str(), "READ" );
  }

  // Exit if no version of file found
  if (infile->IsZombie()){
    cerr << "[ERROR ] : Cannot find flux input file: " << f << endl;
    exit(30);
  } else {
    cout << "Opened flux successfully!" << endl;
  }

  spectrum = (TH1D*) infile->Get(h.c_str());

  // Check histogram okay
  if (!spectrum){
    cerr << "Cannot find histogram : " << h << endl;
    cerr << "File Contents : " << endl;
    infile->ls();
    exit(30);
  }

  // Get Histogram
  spectrum->SetDirectory(0);

  double scaleF = 1.0;
  if (inMEV) scaleF = MeV;
  else scaleF = GeV;

  Emin = spectrum->GetXaxis()->GetXmin() * scaleF;
  Emax = spectrum->GetXaxis()->GetXmax() * scaleF;

  n = spectrum->GetNbinsX();

  prob.resize(n);
  prob2.resize(n);
  Elow.resize(n);
  Ehigh.resize(n);

  double prev1 = 0;
  double prev2 = 0;

  for (int i = 0; i < n; i++){

    double E = spectrum->GetXaxis()->GetBinCenter(i+1) * scaleF;
    Elow[i]  = spectrum->GetXaxis()->GetBinLowEdge(i+1) * scaleF;
    Ehigh[i] = spectrum->GetXaxis()->GetBinLowEdge(i+2) * scaleF;

    prob[i]  = prev1 += spectrum->GetBinContent(i+1);
    prob2[i] = prev2 += spectrum->GetBinContent(i+1) * E;

  }

  // Close Shop
  infile->Close();
}


  
/// yield a random number with the probability profile
/// specified by the constructor parameters
inline double EnergyProfile::shoot(bool dis) {
  if (n == 0) return Emin;
  
  auto &pro = dis ? prob2 : prob;
  double x = frandom() * pro[n - 1];
  int i = 0, j = n - 1;
  while (i < j) {
    int s = (i + j) / 2;
    if (x < pro[s])
      j = s;
    else
      i = s + 1;
  }

  double z = frandom();
  if (dis) {
    double a = Elow[i];
    double b = Ehigh[i];      
    return sqrt(a * a + z * (b * b - a * a));
  }
  return Elow[i] + z * (Ehigh[i] - Elow[i]);
}

/// diagnostic function for verifying if the beam
/// has been correctly constructed
inline void EnergyProfile::print() {
  cout << "Emin=" << Emin << endl;
  cout << "Emax=" << Emax << endl;
  cout << "n=" << n << endl;
  if (n > 1)
    for (int i = 0; i < n; i++)
      cout << "prob[" << i << "]=" << prob[i] - (i ? prob[i - 1] : 0) << endl;
}

inline double EnergyProfile::minE() { return Emin; }

inline TH1D EnergyProfile::GetHist(bool inMEV) const{

  // If not using root
  if (!roothist){
    TH1D temp_flux;

    // Make Histogram
    if (inMEV){
      temp_flux =  TH1D("FluxHist", "FluxHist;E_{#nu} (MeV);Flux",
			n,
			Emin ,
			Emax);
    } else {
      temp_flux=  TH1D("FluxHist", "FluxHist;E_{#nu} (GeV);Flux",
		       n,
		       Emin / 1.E3 ,
		       Emax / 1.E3);
    }

    // Fill Hist
    for (size_t bin = 1; bin < temp_flux.GetXaxis()->GetNbins() + 1; ++bin) {
      temp_flux.SetBinContent(bin, GetBinProb(bin - 1, false));
    }

    return temp_flux;

  } else {

    // Assumes GeV
    TH1D temp_flux = *spectrum;
    temp_flux.SetNameTitle("FluxHist", "FluxHist;E_{#nu} (GeV);Flux");

    return temp_flux;
  }
}
////////////////////////////////////////////////////////////////////////////////
///        END OF  I M P L E M E N T A I O N
////////////////////////////////////////////////////////////////////////////////

#endif
