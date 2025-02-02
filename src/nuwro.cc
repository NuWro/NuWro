#include <iomanip>
#include <sstream>
#include <vector>
#include "event1.h"
#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "qelevent.h"
#include "e_el_event.h"
#include "e_spp_event.h"
#include "hypevent.h"
#include "pdg.h"
#include "chooser.h"
#include "beam.h"
#include "beamHist.h"
#include "target_mixer.h"
#include "stdlib.h"
#include "pauli.h"
#include "dis/singlepion.h"
#include "dis/disevent.h"
#include "dis/resevent2.h"
#include "resevent_hybrid.h"
#include "cohevent2.h"
#include "coh.h"
#include "mecevent.h"
#include "args.h"
#include "kaskada7.h"
#include "sfevent.h"
#include "Analyser1.h"
#include "geomy.h"
#include "ff.h"
#include "hist.h"
#include "nucleusmaker.h"
#include "Interaction.h"
#include "rew/rewparams.h"
#include "lepevent.h"
#include "output.h"


// extern double SPP[2][2][2][3][40];
//extern double sppweight;
extern "C" {
void shhpythiaitokay_(void);
void youcanspeaknowpythia_(void);
}


params *p1=NULL;
string data_dir;
#include "nuwro.h"

NuWro::~NuWro()
{		
	delete _mixer;
	delete _detector;
	delete _beam;
	delete _nucleus;
}

NuWro::NuWro()
{
	_mixer = NULL;
	_detector = NULL;
	_beam = NULL;
	_nucleus = NULL;
}

void NuWro :: set (params &par)
{	
	p = par;
	
	frandom_init(par.random_seed);

	dismode = false;

	if(par.target_type == 1)
		_mixer = new target_mixer (par);
	_detector = make_detector (par);
	
	_beam = create_beam (par,_detector);

	_nucleus = make_nucleus (par);
	
	ff_configure (par);
	refresh_dyn (par);
}

void NuWro :: refresh_target (params &par)
{
	delete _nucleus;
	_nucleus = make_nucleus (par);
}

void NuWro :: refresh_dyn (params &par)
{
	_procesy.reset(par);
}

geomy* NuWro::make_detector(params &p)
{
	if(p.target_type!=2)
		return NULL;
		
	if(p.geo_file.length())
	{
		try
		{
			if(p.geo_d.norm2()==0)
				return new geomy(p.geo_file,p.geo_name,p.geom_length_units,p.geom_density_convert);
			else
				return new geomy(p.geo_file,p.geo_name,p.geom_length_units,p.geom_density_convert,p.geo_volume,p.geo_d,p.geo_o);
		}
		catch(...)
		{
			cerr<<"Failed to make detector."<<endl;
			exit(3);
		}
	}
	else
	{	
		cerr<<"Failed to make detector. Parameter geo_file must not be empty if target_type=2."<<endl;
		exit(4);
		return NULL;
	}

}

void NuWro::init (int argc, char **argv)
{
	frame_top("Simulation parameters");

	//dismode=false;
	dismode=true;
	set_dirs(argv[0]);
	a.read (argc, argv);
	p.read (a.input);
	p.read (a.params, "command line");
	p.list (cout);
	p.list (string(a.output)+".par");
	p1=&p;
	rew.init(p);
	if(a.progress)
		_progress.open(a.progress);
	frandom_init(p.random_seed);

	frame_bottom();

	frame_top("Initialize the simulation");

	if(p.beam_test_only==0 && p.kaskada_redo==0)
		if(p.dyn_dis_nc or p.dyn_res_nc  or p.dyn_dis_cc or p.dyn_res_cc )
	{
		cout << "     -> Calculating the single-pion functions..." << endl;
		singlepion (p);
	}
	if(p.kaskada_redo==0)
	{
		cout << "     -> Building the target nuclei..." << endl;
		_nucleus = make_nucleus(p);
		if(p.target_type==1)
			_mixer=new target_mixer(p);
		else
			_mixer = NULL;
		cout << "     -> Constructing the detector..." << endl;
		_detector=make_detector(p);
        
		cout << "     -> Creating the beam..." << endl;
		_beam=create_beam(p,_detector);
		if(_beam==NULL)
			{
        cerr<<"No beam defined."<<endl;
        exit(5);
      }
    }

  // load the input data
  cout << "     -> Loading external physical data..." << endl;
  input.initialize( p );
  input.load_data();

  cout << "     -> Configuring form factors..." << endl;
  ff_configure(p);

  cout << "     -> Extablishing the choice of dynamics..." << endl;
  refresh_dyn(p);

  frame_bottom();
}

void NuWro::makeevent(event* e, params &p)
{
	static double max_norm=0;
	particle nu;
	int dyn = e->dyn;
	if(_detector)
	{
		material mat;
		do
		{
			nu=_beam->shoot(1<dyn && dyn<6 && dismode);
			if(nu.travelled>0 && p.beam_weighted==0)
			{
				if(nu.travelled<frandom()*max_norm)
					continue;
				if(nu.travelled>max_norm)
					max_norm=nu.travelled;
			}			
			nu.travelled=0;
			nu.r=vec(nu.r)+p.beam_offset;
			if(nu.r.x==0 && nu.r.y==0 && nu.r.z==0)
				mat=_detector->getpoint();
			else
				mat=_detector->getpoint(nu.p(),nu.r);

		} while(not (mat.Z+mat.N>0 && mat.w_density>=frandom()*_detector->max_dens()));

		///change nucleus
		e->r=mat.r;
		p.nucleus_p=mat.Z;
		p.nucleus_n=mat.N;
		p.nucleus_E_b=0;  // Use library value for E_b
		p.nucleus_kf=0;   // Use library value for kF
//		cout<<mat.Z<<' '<<mat.N<<' '<<endl;
		if(mat.Z==0&&mat.N==0)
			throw "Empty isotope 00";
	}
	else
	{
		nu=_beam->shoot(1<dyn && dyn<6 && dismode);
		nu.r=vec(nu.r)+p.beam_offset;
	}

	if(_detector or _mixer) // _nucleus not reusable
	{
		delete _nucleus;
		_nucleus= make_nucleus(p);
		//cout<<"make_nucleus "<<_nucleus->p<<" "<<_nucleus->n<<endl;
	}
	else
		_nucleus->reset();
	e->in.push_back (nu);		 // insert neutrino
	if(dyn<6 || (dyn>=10 && dyn<12) || dyn==20)
	{
								 // insert target nucleon
		e->in.push_back (_nucleus->get_nucleon());
		e->in[0].r=e->in[1].r;
		assert(e->in[1]*e->in[1]>0);
	}
	else if(dyn>=12 && dyn<14)
	{
		// insert target electron
		e->in.push_back(particle(PDG::pdg_e, PDG::mass_e));
	}

	e->weight=0;
	if(nu.travelled>0)
		e->norm=nu.travelled;
	// else e->norm remains 1;

	e->flag.cc  = false;
	e->flag.nc  = false;
				  
	e->flag.qel = false;
	e->flag.res = false;
	e->flag.dis = false;
	e->flag.coh = false;
	e->flag.mec = false;
	e->flag.hyp = false;
	e->flag.lep = false;
	
	e->flag.anty = nu.pdg<0;

	if(p.beam_test_only)
	{
		e->weight=1;
		e->out.push_back(e->in[0]);
		return;
	}
	double factor=1.0;
								 
	if(p.cc_smoothing and dyn==0) //only in qel_cc
	{
		if(e->in[0].pdg>0)
		{
			factor=_nucleus->frac_neutron();
			e->in[1].set_neutron();
		}
		else
		{
			factor=_nucleus->frac_proton();
			e->in[1].set_proton();
		}
	}
	e->par =p;
	
	
	if(  // (anty)-neutrino interaction
	   abs(nu.pdg)==12 or abs(nu.pdg)==14 or abs(nu.pdg)==16
	  )
	switch (dyn)
	{
		case 0: 
			e->flag.qel=e->flag.cc=true;
			if (p.dyn_qel_cc) // qel cc
			{
				if(p.sf_method>0 and has_sf(*_nucleus, p.sf_method))
					sfevent (p, *e, *_nucleus);
				else
					qelevent1 (p, *e, *_nucleus, false);
			}
			break;				 
		case 1:
			e->flag.qel=e->flag.nc=true;
			if (p.dyn_qel_nc) // qel nc
			{
				if(p.sf_method>0 and has_sf(*_nucleus, p.sf_method))
					sfevent (p, *e, *_nucleus);
				else
				qelevent1 (p, *e, *_nucleus, true);
			}
			break;				 
		case 2:
			e->flag.res=e->flag.cc=true;
			if (p.dyn_res_cc) // res cc
			{
				switch(p.res_kind)
				{
					case 1:resevent2 (p, *e, *_nucleus, true);break;
					case 2:resevent_hybrid (p, *e, *_nucleus, true);break;
					default:resevent2 (p, *e, *_nucleus, true);break;
				}
				if (p.pauli_blocking)
					mypauli_spp (*e, *_nucleus);
			}
			break;				
		case 3:
			e->flag.res=e->flag.nc=true;
			if (p.dyn_res_nc) // res nc
			{
				resevent2 (p, *e, *_nucleus, false);
				if (p.pauli_blocking)
					mypauli_spp (*e, *_nucleus);
			}
			break;
		case 4:
			e->flag.dis=e->flag.cc=true;
			if (p.dyn_dis_cc) // dis cc
			{
				disevent (p, *e, *_nucleus, true);
				if (p.pauli_blocking)
					mypauli_spp (*e, *_nucleus);
			}
			break;				
		case 5:
			e->flag.dis=e->flag.nc=true;
			if (p.dyn_dis_nc) //dis nc
			{
				disevent (p, *e, *_nucleus, false);
				if (p.pauli_blocking)
					mypauli_spp (*e, *_nucleus);
			}
			break;
		case 6:                  
			e->flag.coh=e->flag.cc=true;
			if (p.dyn_coh_cc) // coh cc
			{
				if(p.coh_new)
				switch(p.coh_kind)
				{
					case 1:cohevent_cj (p, *e, *_nucleus, true);break;
					case 2:cohevent_bs (p, *e, *_nucleus, true);break;
					default:cohevent_bs (p, *e, *_nucleus, true);break;
				}
				else          cohevent2   (p, *e, *_nucleus, true);
			}
			break;
		case 7:                  
			e->flag.coh=e->flag.nc=true;
			if (p.dyn_coh_nc) // coh nc
			{
				if(p.coh_new)
				switch(p.coh_kind)
				{
					case 1:cohevent_cj (p, *e, *_nucleus, false);break;
					case 2:cohevent_bs (p, *e, *_nucleus, false);break;
					default:cohevent_bs (p, *e, *_nucleus, false);break;
				}
				else          cohevent2   (p, *e, *_nucleus, false);
			}
			break;
		case 8:
			e->flag.mec=e->flag.cc=true;
			if (p.dyn_mec_cc) // mec cc
			//if( nu.pdg>0 || !(p.mec_kind==3) )// al flavor states/antineutrinos available
			{
				if(_nucleus->A()<=1)
					break;
				switch(p.mec_kind)
				{
					case 1:mecevent_tem (p, *e, *_nucleus, true);break;
					case 2:mecevent2 (p, *e, *_nucleus, true, false);break;
					case 3:mecevent_Nieves (p, *e, *_nucleus, true);break;
					case 4:mecevent2 (p, *e, *_nucleus, true, true);break;
					case 5:mecevent_SuSA (p, *e, *_nucleus, true);break;
					default:mecevent_tem (p, *e, *_nucleus, true);break;
				}
				for(int i=0;i<e->out.size();i++)
				{
					e->out[i].r=e->in[1].r;
					e->out[i].set_momentum(e->out[i].p().fromZto(e->in[0].p()));
				}
			}
			break;
		case 9:
			e->flag.mec=e->flag.nc=true;
			if (p.dyn_mec_nc) //mec nc
			if(p.mec_kind==1)      // only TEM for NC
			{
				if(_nucleus->A()<=1)
					break;
				switch(p.mec_kind)
				{
					case 1: mecevent_tem(p, *e, *_nucleus, false);break;
					default: mecevent_tem (p, *e, *_nucleus, false);break; 
				}
				for(int i=0;i<e->out.size();i++)
				{
					e->out[i].r=e->in[1].r;
					e->out[i].set_momentum(e->out[i].p().fromZto(e->in[0].p()));
				}
			}
			break;
		case 10:
			e->flag.hyp=e->flag.cc=true;
			if(p.dyn_hyp_cc) // qel hyperon
			{
				hypevent (p, *e, *_nucleus);
			}
			break;		
		case 12:
			e->flag.lep=true; //->flag.cc=true;
			if (p.dyn_lep) // Neutrino-lepton
			{
				lepevent (p, *e); //, true);
			}
			break;
	}
	else if(e->in[0].pdg==11) // electron scattering
	{
       	switch(dyn)
		{
			case 20:
				// TODO: introduce a new flag el!
				e->flag.qel=e->flag.nc=true;
			    /*if(p.eel_alg=="old")
                    e_el_event(p,*e,*_nucleus,false); 
			    else 	
                if(p.eel_alg=="fast")
                    e_el_event2orig(p,*e,*_nucleus,false); 
                else   // all remaining algorithms  
                    e_el_event2(p,*e,*_nucleus,false); */
                if(p.sf_method>0 and has_sf(*_nucleus, p.sf_method))
					sfevent (p, *e, *_nucleus);
				else
					qelevent1 (p, *e, *_nucleus, true);
			    break;
			// case 21: 
			// 	e->flag.nc=true;
			//     if(p.eel_theta_lab>0) 	
   //                  e_spp_event(p,*e,*_nucleus,false); 
			//     else // use negative theta to test new implementation	
   //                  e_spp_event3(p,*e,*_nucleus,false); 
			//     break;
		}
	}
	e->weight*=factor;

	if (e->weight == 0)
	{
		e->out.clear ();
								
		e->out.push_back (e->in[0]);    
		e->out.push_back (e->in[1]);
	}
	//      e->check();
}	// end of makeevent


void NuWro::finishevent(event* e, params &p)
{
  // Resample independent variables not used for event acceptance
  if( e->flag.res && e->flag.cc )
  {
    if( p.res_kind == 2 && e->flag.res_delta )
    {
      // for consistency reasons, the simplest solution for now is to recreate the nucleus
      nucleus *nucl = make_nucleus(p);
      if( p.res_hybrid_sampling == 1 && e->flag.need_resample_dir )
        resevent_dir_hybrid(*e, *nucl, p.res_hybrid_resampling);
      else if( p.res_hybrid_sampling < 4 && e->flag.need_resample_phi )
        resevent_phi_hybrid(*e, *nucl);
      delete nucl;
    }
  }

	for(int i=0;i<1/* e->in.size()*/;i++)
	{
		e->in[i].endproc=e->dyn;
		registration(e->all,e->in[i]);
	}
	for(int i=0;i<e->out.size();i++)
	{
		e->out[i].mother=0;
		registration(e->all,e->out[i]);
	}
	if(p.beam_test_only)
		return;

	for(int j=0;j<e->in.size();j++)
	{
		particle p=e->in[j];
		p.endproc=e->dyn;
		registration(e->all,p);
	}

	//e->pr=_nucleus->Zr(); 	// 1. po co to?
	//e->nr=_nucleus->Nr(); 	// 2. powoduje break, segmentation fault

								 // copy particle from out to post if coherent interaction
	
	if (!e->flag.coh && !e->flag.lep && (e->par.nucleus_n + e->par.nucleus_p > 1))
	{
		kaskada k(p, *e, &input);
		k.kaskadaevent();		 // runs only if p.kaskada_on is true
	}
	else
//	if(e->post.size()==0)   // copy out to post if no fsi
	{
		for (int j = 0; j<e->out.size(); j++)
		{
			particle p = e->out[j];
			p.endproc = e->dyn;
			registration(e->all,p);
			e->post.push_back(p);
		}
	}
}								 //end of finishevent

void NuWro::raport(double i, double n, const char* text, int precision, int k, string label, bool toFile)
{
  static int prev=-1;
  int proc=precision*i/n;
  if(proc!=prev)
  {
    prev=proc;
    cerr.precision(3);
    if(toFile)
    {
      _progress.seekp(ios_base::beg);
      _progress << proc*100.0/precision << " " << text << '\r' << flush;
      _progress.flush();
    }
    else
    {
      cerr << "        ";
      cerr << showpoint << proc*100.0/precision << " % of ";
      if(k>=0)
        cerr << label << " ";
      cerr << text << '\r' << flush;
    }
    cerr.precision();
  }
}

void NuWro::pot_report(ostream& o, bool format=false)
{
	double tot=0;
	if(_detector)
	{
		if(format)
			frame_top("POT Report");
		string tab(8,' ');
		if(!format)
			tab.clear();
		double pd=_detector->nucleons_per_cm2();

		for(int i=0;i<_procesy.size();i++)
		if(_procesy.avg(i)>0)
		{	
			tot+=_procesy.avg(i);
			double epp=_procesy.avg(i)*pd*_beam->nu_per_POT();
			o  <<tab<<"dyn["<<i<<"] events per POT = "<<epp<<endl;
			o  <<tab<<"       POT per event = "<<1/epp<<endl;
		}
		double epp=tot*pd*_beam->nu_per_POT();
	
		o<<tab<<"Total: events per POT= "<<epp<<endl
		 <<tab<<"       POT per event = "<<1.0/epp<<endl;
		o<<tab<<" BOX nuclons per cm2 = "<<pd<<endl;
		o<<tab<<"Total cross section  = "<<tot<<" cm2"<<endl;
		o<<tab<<"Reaction probability = "<<tot*pd<<endl;
		o<<tab<<"Average BOX density  = "<<_detector->density()/g*cm3<<" g/cm3"<< endl;		
		o<<tab<<"Estimated BOX mass   = "<<_detector->vol_mass()/kg<<" kg"<<endl;	
		o<<tab<<"Fraction of protons  = "<<_detector->frac_proton()<<endl;
		o<<tab<<"Total POT: " << p.number_of_events/epp << endl << endl;
		if(format)
			frame_bottom();
	}
}


//////////////////////////////////////////////////////////////
//              Test events
//////////////////////////////////////////////////////////////
void NuWro::test_events(params & p)
{
	frame_top("Run test events");

	if(p.number_of_test_events>0  && p.beam_test_only==0)
	{
		hist hq2((char*)"q2",0,2*GeV2,100);
		hist hq0((char*)"q0",0,(p.beam_type==0 ? atof(p.beam_energy.c_str())*MeV : 2*GeV),100);
		hist hqv((char*)"qv",0,(p.beam_type==0 ? atof(p.beam_energy.c_str())*MeV : 2*GeV),100);
		hist hT((char*)"T",0,2*GeV,100);
		TFile *te=NULL;
		TTree *t1=NULL;
		event *e=NULL;
		if(p.save_test_events)
		{	
			dismode=false;
			te=new TFile((string("weighted.")+a.output).c_str(),"recreate");						
			t1 = new TTree ("treeout", "Tree of events");
			e = new event ();
			t1->Branch ("e", "event", &e);
			delete e;
		}

		refresh_dyn(p);
		  
		int saved=0;
		for (int i = 0; i < p.number_of_test_events; i++)
		{
			e = new event ();
			int k= _procesy.choose(); 
			e->dyn = _procesy.dyn(k); // choose dynamics
			if(_mixer)
				_mixer->prepare(p);
			makeevent(e,p);
			double bias=1;
			if(dismode && e->dyn>1 && e->dyn<6)
				bias=e->in[0].t;

			_procesy.add (k, e->weight, bias);
			e->weight/=_procesy.ratio(k); // make avg(weight)= total cross section
			//~ if(_detector and _beam->nu_per_POT() != 0)
			//~   e->POT=e->weight * _detector->nucleons_per_cm2() / _beam->nu_per_POT();

			if(e->weight>0)
			{
				hq2.insert_value(-e->q2(),e->weight*cm2);
				hq0.insert_value(e->q0(),e->weight*cm2);
				hqv.insert_value(e->qv(),e->weight*cm2);
				hT.insert_value(e->in[0].E(),e->weight*cm2);
			}
			else
			{
				hq2.insert_value(0,0);
				hq0.insert_value(0,0);
				hqv.insert_value(0,0);
				hT.insert_value(0,0);
			}
			switch(p.save_test_events)
			{
				case 0: 
					break;
				case 1: 
					// finishevent(e, p);
					t1->Fill ();
					break;
				case 2:
					if(e->weight>0)
					{
						saved++;
						e->weight=e->weight*saved/(i+1);
						finishevent(e, p);
						t1->Fill ();						
					}
					break;
				default:
					cerr<<"Parameter save_test_events="<<p.save_test_events;
					cerr <<" out of range. Should be ,1, or 2)"<<endl;
					exit(12);
			}
			delete e;
			raport(i+1,p.number_of_test_events,"test events ready... ",1000,-1,"",bool(a.progress));
		}						 // end of nuwro loop
		if(p.save_test_events)
		{
			te->Write ();
			te->Close ();
		}

    cout << "        100. % of test events ready..." << endl;

		_procesy.report();
		_procesy.set_weights_to_avg ();
		string prefix;
//		if(strlen(a.output)>5 && string(".root")==a.output[strlen(a.output)-5])
			prefix="";
//		else 
//			prefix=a.output;
		hq2.plot(prefix+"q2.txt",GeV2,1e-38*cm2/GeV2);
		hq0.plot(prefix+"q0.txt",GeV,1e-38*cm2/GeV);
		hqv.plot(prefix+"qv.txt",GeV,1e-38*cm2/GeV);
		hT.plot(prefix+"T.txt",GeV,1e-38*cm2/GeV);
		
		ofstream totals ((prefix+"totals.txt").c_str(),ios::app);
		totals<<p.beam_energy;
		double tot=0;
		
		int j=0;
		for(int i=0;i<_procesy.size();i++)
		{
			while(j++<_procesy.dyn(i))
				totals << ' '<<0; // cross section of disabled channels
			totals << ' '<<_procesy.avg(i);
			tot+=_procesy.avg(i);
		}
		while(j++<10)
			totals << ' '<<0; // cross section of disabled channels

		totals<<endl;
		pot_report(cout,true);		
		if(_detector)
		{   ofstream  potinfo("POTinfo.txt");
			pot_report(potinfo);
		}
	}
}


void NuWro::user_events(params &p)
{
	if(p.number_of_test_events<1 or p.user_events==0)
		return;
	frame_top("Run analyser events");
	params p1=p;
	Analyser * A=make_analyser(p);
	if(!A) 
		return;
	for(A->start(); !A->end(); A->step())
	{
		refresh_dyn(p);
		for (int i = 0; i < p.number_of_test_events; i++)
		{
			event *e = new event ();
								 
			int k = _procesy.choose ();///< choose bin
			e->dyn = _procesy.dyn(k);  ///< choose dynamics
			
			A->prepare_event(*e);
			if(_mixer)
				_mixer->prepare(p);

			makeevent(e,p);

			A->process_event(*e);

			double bias=1;
			if(dismode && e->dyn>1 && e->dyn<6)
				bias=e->in[0].t;

			_procesy.add (k, e->weight, bias);

			delete e;
			
			raport(i+1,p.number_of_test_events,"analyser events ready... ",1000,-1,"",bool(a.progress));
			p=p1;
	}	// end of nuwro loop

    cout << "        100. % of analyser events ready..." << endl;

		A->partial_report();
		_procesy.report();
	
	}	// end of analyser loop
	A->final_report();
	delete A;
}


void NuWro::real_events(params& p)
{
	dismode=true;
	if(p.number_of_events<1)
		return;

	frame_top("Update channels ratio");

	/// calculate desired number of events for each dynamics
	_procesy.calculate_counts(p.number_of_events);
	{							 /// Write cross sections and counts to screen and file
		ofstream f((string(a.output)+".txt").c_str());
		_procesy.short_report(cout,true);
		_procesy.short_report(f);
	}

  	frame_bottom();

	frame_top("Run real events");

	event *e = new event;

	string output=a.output;
	int l=output.length();
	if(l<5 || string(".root")!=output.c_str()+l-5)
		output=output+".root";
	TFile *ff = new TFile (output.c_str(), "recreate");
	TTree *tf = new TTree ("treeout", "Tree of events");
	tf->Branch ("e", "event", &e);
	delete e;
	TH1 * xsections= new TH1D("xsections","xsections",_procesy.size(),0,_procesy.size());
	for(int i=0;i<_procesy.size();i++)
	{
		xsections->SetBinContent(i+1,_procesy.avg(i));
		xsections->SetBinError(i+1,_procesy.sigma(i));
	}

	//	TNamed version("NuWro version", VERSION );
	//	version.Write();

	/////////////////////////////////////////////////////////////
	// The main loop in NPROC -- generate file with unweighted events
	//////////////////////////////////////////////////////////////

	char filename[230];
	if(p.kaskada_redo==0)
		for (int k = 0; k < _procesy.size(); k++)
			if(_procesy.desired(k)>0)
			{
				sprintf(filename,"%s.%d.part",a.output,k);
				TFile *f1 = new TFile (filename, "recreate");
				TTree *t1 = new TTree ("treeout", "Tree of events");

				e = new event ();
				t1->Branch ("e", "event", &e);
				delete e;

				while(_procesy.ready(k)<_procesy.desired(k))
				{
					e = new event ();
					e->dyn = _procesy.dyn(k);

					if(_mixer)
						_mixer->prepare(p);
					makeevent(e,p);
					double bias=1;
					if(!p.beam_test_only && dismode && e->dyn>1 && e->dyn<6)
						bias=e->in[0].t;
					if (_procesy.accept(k,e->weight,bias))
					{
						finishevent(e, p);
						e->weight=_procesy.total();
						//~ if(_detector and _beam->nu_per_POT() != 0)
							//~ e->POT=e->weight * _detector->nucleons_per_cm2() / _beam->nu_per_POT();
						t1->Fill ();
					}
					delete e;

					raport(_procesy.ready(k),_procesy.desired(k),"events ready... ",1000,_procesy.dyn(k),_procesy.label(k),bool(a.progress));
				}
				f1->Write ();

        cout << "        100. % of " << _procesy.label(k) << " events ready..." << endl;

		// elimination of spurious events for dynamics k
		// by copying only last desired[k] events to outfile
				if(p.mixed_order==0)
				{
					int nn = t1->GetEntries ();
					int start = nn-_procesy.desired(k);
					for (int jj = start; jj < nn; jj++)
					{
						e = new event();
						t1->GetEntry (jj);
						tf->Fill ();
						delete e;
						raport(jj-start+1,nn-start,"events copied... ",100,_procesy.dyn(k),_procesy.label(k),bool(a.progress));
					}
				}

				f1->Close ();
				delete f1;
				if(p.mixed_order==0)
					unlink(filename);
			};
	//////////////////////////////////////////////////////////////////////////////////////
	//                    end of the main loop in NPROC
	//////////////////////////////////////////////////////////////////////////////////////
	if(p.mixed_order)
	{
		TFile *f[_procesy.size()];
		TTree *t[_procesy.size()];
		int n[_procesy.size()],u[_procesy.size()];
		int ile=0;
		e=new event();
		for (int k = 0; k < _procesy.size(); k++)
			if((u[k]=_procesy.desired(k))>0)
		{
			sprintf(filename,"%s.%d.part",a.output,k);
			f[k] = new TFile (filename);
			t[k] = (TTree*) f[k]->Get("treeout");
			if(t[k]==NULL)
			{
				cerr<< "tree \"treeout\" not found in file \""<<filename<<"\""<<endl;
				exit(6);
			}

			t[k]->SetBranchAddress("e", &e);
			n[k]=t[k]->GetEntries();
			ile+=u[k];
		}
		int nn=ile;
		while(ile>0)
		{
			int i=0;
			int x=frandom()*ile;
			while(x>=u[i])
				x-=u[i++];
			t[i]->GetEntry(n[i]-u[i]);
			tf->Fill();
			ile--;
			u[i]--;
			raport(nn-ile,nn,"real events copied... ",100,-1,"",bool(a.progress));
		}
		delete e;
		for (int k = 0; k < _procesy.size(); k++)
			if(_procesy.desired(k))
		{
			f[k]->Close();
			delete f[k];
			sprintf(filename,"%s.%d.part",a.output,k);
			unlink (filename);
		}
	}
	ff->Write ();
	ff->Close ();
	cout << "        100. % of real events copied..." << endl;
	_progress.close();
	delete ff;
	_procesy.report();
	pot_report(cout,true);
	if(_detector)
	{   ofstream  potinfo("POTinfo.txt");
		pot_report(potinfo);
	}

	frame_top("Finalize the simulation");
	cout << "     " << "-> Generated the output file: \"" << output << "\"" << endl;
	frame_bottom();
}


void NuWro::kaskada_redo(string input,string output)
{
	event *e=new event;

	TFile *fi = new TFile (input.c_str());
	TTree* ti = (TTree*) fi->Get("treeout");
	if(ti==NULL)
	{
		cerr<< "tree \"treeout\" not found in file \""<<input<<"\""<<endl;
		exit(7);
	}
	ti->SetBranchAddress("e", &e);

	TFile *ff= new TFile(output.c_str(),"recreate");
	TTree *tf = new TTree("treeout","Tree of events");
	tf->Branch("e","event",&e);

	int nn = ti->GetEntries ();
	for (int i = 0; i < nn; i++)
	{
		//		e = new event();
		ti->GetEntry (i);
		e->clear_fsi();
		finishevent(e,p);
		tf->Fill ();
		//		delete e;
		//if(i%1000==0)
		//cout<<i/1000<<"/"<<nn/1000<<"\r"<<endl;
		raport(i+1,nn," % events processed...",100,e->dyn,_procesy.label(e->dyn),bool(a.progress));
	}
	cout<<endl;
	fi->Close ();
	delete fi;

	ff->Write();
	ff->Close();
	delete e;
	delete ff;
	cout<<"Output: \""<<output<<"\""<<endl;
}


#include "UserAction.h"

void NuWro::main (int argc, char **argv)
{
  shhpythiaitokay_();
	try
	{
		init(argc,argv);
		if(p.kaskada_redo==1)
			kaskada_redo(a.output,string(a.output)+".fsi.root");
		else
		{
			if(not p.beam_test_only)
			{
				if(p.user_events>0)
					user_events(p);
				else
				{
					test_events(p);
					real_events(p);
				}
			}
		}
		genrand_write_state();
	}
	catch(string s)
	{
		cout<<s<<endl;
	}
	catch(char const* s)
	{
		cout<<s<<endl;
	}
	catch(...)
	{
		cout<<"NuWro failed"<<endl;
	}
}
