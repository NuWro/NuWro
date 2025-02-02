#include "mb_nce.h"
#include "TMinuit.h"

const bool waga = false;

void sim()
{	
	string resF = "mb_nce/root/C_res.root";
	string qelH = "mb_nce/root/H_qel.root";
	
	if (noFile(resF))
	{
		string command = "./bin/nuwro -o '" + resF +"' " + parRES + fg + ccma + ncma(1350) + C + wf(7);
		run(command);
	}
	
	if (noFile(qelH))
	{	
		string command = "./bin/nuwro -o '" + qelH + "' " + parQEL + fg + ccma + ncma(1030) + H;
		run(command);
	}
	
	for (int i = 0; i < nof_ma; i++)
	{				
		stringstream temp;
		string x;
	
		temp << ma[i];
		temp >> x;
		
		string qelC = "mb_nce/root/C_qel_" + x + ".root";
		
		if (noFile(qelC))
		{
			string command = "./bin/nuwro -o '" + qelC + "' " + parQEL + fg + ccma + ncma(ma[i]) + C + wf(7);
			run(command);
		}
	}
}

void re_sim()
{
	string resF = "mb_nce/root/C_res.root";
	
	string command = "./bin/nuwro -o '" + resF + "' -p 'kaskada_redo = 1' ";
	run(command);
		
	run("mv " + resF + ".fsi.root " + resF);
		
	for (int i = 0; i < nof_ma; i++)
	{			
		stringstream temp;
		string x;
		
		temp << ma[i];
		temp >> x;
	
		string qelC = "mb_nce/root/C_qel_" + x + ".root";
	
		string command = "./bin/nuwro -o '" + qelC + " -p 'kaskada_redo = 1' ";
		run(command);
			
		run("mv " + qelC + ".fsi.root " + qelC);
	}
}

void calc()
{	
	string qelH = "mb_nce/root/H_qel.root";
	string qelHtxt = qelH + ".txt";
	
	double xsecH = crosssection(qelHtxt);
	
	if (waga)
		xsecH = 1.0;
	
	double Hresult[51] = {0};
	calcH(qelH, Hresult);
	
	string rdC = "mb_nce/root/C_res.root";
	string rdCtxt = rdC + ".txt";

	double result[5][51] = {{0}};
		
	calcC(rdC, result);
	double resxsec = crosssection(rdCtxt);
	
	if (waga)
		resxsec = 1.0;

	for (int i = 0; i < 51; i++)
	{
		result[4][i] *= resxsec/events/18.0;
		result[0][i] = Hresult[i]*xsecH/events/18.0;
	}
		
	norm(result[0], 2.0/14);
	norm(result[4], 12.0/14);
		
	for (int j = 0; j < nof_ma; j++)
	{
		for (int i = 0; i < 51; i++)
		{
			result[1][i] = 0;
			result[2][i] = 0;
			result[3][i] = 0;
		}

		stringstream temp;
		string x;
		
		temp << ma[j];
		temp >> x;

		string qelC = "mb_nce/root/C_qel_" + x + ".root";
		string qelCtxt = qelC + ".txt";
		
		string out = "mb_nce/results/true_" + x + ".txt";
		string out2 = "mb_nce/results/rec_" + x + ".txt";
		
		if (!noFile(out) and !noFile(out2))
			continue;
			
		calcC(qelC, result);
			
		double xsecC = crosssection(qelCtxt);
		
		if (waga)
			xsecC = 1.0;
			
		for (int i = 0; i < 51; i++)
		{
			result[1][i] *= xsecC/events/18.0;
			result[2][i] *= xsecC/events/18.0;
			result[3][i] *= xsecC/events/18.0;
		}
			
		norm(result[1], 12.0/14);
		norm(result[2], 12.0/14);
		norm(result[3], 12.0/14);
			
		double result2[51];
			
		true2rec(result, result2);
			
		ofstream file(out.c_str());
			
		for (int i = 0; i < 51; i++)
		{
			file << i*18 << " " << result[0][i] << " " << result[1][i] << " " << result[2][i] << " " << result[3][i] << " " << result[4][i] << endl;
			file << (i+1)*18 << " " << result[0][i] << " " << result[1][i] << " " << result[2][i] << " " << result[3][i] << " " << result[4][i] << endl;
		}
			
		file.close();

		ofstream file2(out2.c_str());
			
		for (int i = 0; i < 51; i++)
		{
			file2 << 40 + i*12 << " " << result2[i] << endl;
			file2 << 40 + (i+1)*12 << " " << result2[i] << endl;
		}
		
		file2.close();
	}
}

void calcH(string in, double *res)
{
	TFile *tf1 = new TFile(in.c_str());
	TTree *tt1 = (TTree*)tf1->Get("treeout");
	event *e1   = new event();
		
	tt1->SetBranchAddress("e",&e1);
	
	for (int i = 0; i < events; i++)
	{
		tt1->GetEntry(i);
		
		if (e1->flag.qel == 1)
		{		
			double Tk = e1->out[1].Ek();
			int w = Tk/18;
			if (w > 50) w = 50;
			
			if (waga)
				res[w] += e1->weight;				
			else
				res[w]++;
		}
		
		cout << in << ": " << 100*i/events << "%\r" << flush;
	}
	
	cout << in << ": done" << endl << endl;
		
	delete e1;
	delete tt1;
	delete tf1;
}

void calcC(string in, double res[5][51])
{
	TFile *tf1 = new TFile(in.c_str());
	TTree *tt1 = (TTree*)tf1->Get("treeout");
	event *e1   = new event();
		
	tt1->SetBranchAddress("e",&e1);
	
	for (int i = 0; i < events; i++)
	{		
		tt1->GetEntry(i);
		
		int piony = e1->fof(211) + e1->fof(-211) + e1->fof(111);
		int nukleony = e1->fof(2212) + e1->fof(2112);
		
		if (piony == 0 and nukleony != 0)
		{			
			double Tk = 0;
			
			for (int k = 0; k < e1->post.size(); k++)
				if (e1->post[k].pdg == 2112 or e1->post[k].pdg == 2212)
					Tk += e1->post[k].Ek();
						
			int a = Tk/18;
						
			if (a > 50)
				a = 50;
			else if (a < 0)
				continue;

			int b = 4;
			
			if (e1->flag.qel)
			{
				if (e1->in[1].pdg == 2112)
					b = 3;
				else if (e1->number_of_interactions() == 0)
					b = 1;
				else
					b = 2;
			}

			if (waga)
				res[b][a] += e1->weight;
			else
				res[b][a]++;
		}
		
		cout << in << ": " << 100*i/events << "%\r" << flush;
	}
	
	cout << in << ": done" << endl << endl;
		
	delete e1;
	delete tt1;
	delete tf1;
}

void delta_s_sim(double xMa)
{	
	for (int i = 0; i < nof_ds; i++)
	{
		stringstream temp;
		string x;
	
		temp << xds[i];
		temp >> x;
		
		string resC = "mb_nce/root/C_res_s" + x + ".root";
		string qelH = "mb_nce/root/H_qel_s" + x + ".root";
		string qelC = "mb_nce/root/C_qel_s" + x + ".root";
		
		if (noFile(resC))
		{
			string command = "./bin/nuwro -o '" + resC + "' " + parRES + fg + ccma + ncma(1350) + C + wf(7) + ds(xds[i]);
			run(command);
		}

		if (noFile(qelH))
		{
			string command = "./bin/nuwro -o '" + qelH + "' " + parQEL + fg + ccma + ncma(1030) + H + ds(xds[i]);
			run(command);
		}
		
		if (noFile(qelC))
		{
			string command = "./bin/nuwro -o '" + qelC + "' " + parQEL + fg + ccma + ncma(xMa) + C + wf(7) + ds(xds[i]);
			run(command);
		}
	}
}

void delta_s_calc()
{	
	for (int i = 0; i < nof_ds; i++)
	{
		stringstream temp;
		string x;
	
		temp << xds[i];
		temp >> x;
		
		string Hin = "mb_nce/root/H_qel_s" + x + ".root";
		string Cin = "mb_nce/root/C_qel_s" + x + ".root";
		string resin = "mb_nce/root/C_res_s" + x + ".root";
		
		string Htxt = Hin + ".txt";
		string Ctxt = Cin + ".txt";
		string restxt = resin + ".txt";
		
		double xsecH = crosssection(Htxt);
		double xsecC = crosssection(Ctxt);
		double xsecres = crosssection(restxt);
		
		if (waga)
		{
			xsecH = 1.0;
			xsecC = 1.0;
			xsecres = 1.0;
		}
			
		double res[5][30] = {{0}};
		
		calcH_he(Hin, res[0]);
		calcC_he(resin, res);
		calcC_he(Cin, res);
		
		for (int j = 0; j < 5; j++)
		{
			double factor = 12.0/14.0*xsecC;
			
			if (j == 0)
				factor = 2.0/14.0*xsecH;
			else if (j == 4)
				factor = 12.0/14.0*xsecres;
			
			normHE(res[j], factor/events/600.0*28.0);
		}
		
		cout << endl;
		
		double singlep[30] = {0};
		double multip[30] = {0};
		
		true2recp(res, singlep);
		true2recN(res, multip);
		
		string out = "mb_nce/results/single_multi_ratio_s" + x + ".txt";
		string out2 = "mb_nce/results/single_multi_ratio_s" + x + "_true.txt";
		
		ofstream file(out.c_str());
		ofstream file2(out2.c_str());
		
		for (int i = 0; i < 30; i++)
		{
			file << 350 + i*15 << " " << singlep[i]/multip[i] << " " << singlep[i] << " " << multip[i] << endl;
			file << 350 + (i+1)*15 << " " << singlep[i]/multip[i] << singlep[i] << " " << multip[i] << endl;
			
			file2 << 300 + i*600.0/28 << " " << res[0][i] << " " << res[1][i] << " " << res[2][i] << " " << res[3][i] << " " << res[4][i] << endl;
			file2 << 300 + (i+1)*600.0/28 << " " << res[0][i] << " " << res[1][i] << " " << res[2][i] << " " << res[3][i] << " " << res[4][i] << endl;
		}
		file.close();
		file2.close();
	}
}

void calcH_he(string in, double *res)
{
	TFile *tf1 = new TFile(in.c_str());
	TTree *tt1 = (TTree*)tf1->Get("treeout");
	event *e1   = new event();
		
	tt1->SetBranchAddress("e",&e1);
	
	for (int i = 0; i < events; i++)
	{
		tt1->GetEntry(i);
		
		if (e1->flag.qel == 1)
		{	
			double cos = e1->out[1].p().z/e1->out[1].momentum();
			
			if (cos <= 0.5)
				continue;
						
			double Tk = e1->out[1].Ek();
			double bin = 600.0/28.0;
			int w = Tk/bin - 13;
			if (w > 29) w = 29;
			else if (w < 0) w = 0;
			
			if (waga)
				res[w] += e1->weight;			
			else
				res[w]++;
		}
		
		cout << in << ": " << 100*i/events << "%\r" << flush;
	}
	
	cout << in << ": done" << endl << endl;
		
	delete e1;
	delete tt1;
	delete tf1;
}

void calcC_he(string in, double res[5][30])
{
	TFile *tf1 = new TFile(in.c_str());
	TTree *tt1 = (TTree*)tf1->Get("treeout");
	event *e1   = new event();
		
	tt1->SetBranchAddress("e",&e1);
	
	for (int i = 0; i < events; i++)
	{		
		tt1->GetEntry(i);
		
		int piony = e1->fof(211) + e1->fof(-211) + e1->fof(111);
		int nukleony = e1->fof(2212) + e1->fof(2112);
		
		if (piony == 0 and nukleony != 0)
		{			
			double Tk = 0;
			
			for (int k = 0; k < e1->post.size(); k++)
				if (e1->post[k].pdg == 2112 or e1->post[k].pdg == 2212)
					Tk += e1->post[k].Ek();
			
			double bin = 600.0/28.0;
			int a = Tk/bin - 13;
			
			if (a > 29) a = 29;			
			else if (a < 0) a = 0;

			int b = 4;
			
			if (e1->flag.qel)
			{
				if (e1->in[1].pdg == 2112)
					b = 3;
				else if (e1->number_of_interactions() == 0)
					b = 1;
				else
					b = 2;
			}
			
			double cos = 1;
			
			if (b == 1)
				cos = e1->out[1].p().z/e1->out[1].momentum();

			if (cos > 0.5)
			{
				if (waga)
					res[b][a] += e1->weight;
				else
					res[b][a]++;
			}
		}
		
		cout << in << ": " << 100*i/events << "%\r" << flush;
	}
	
	cout << in << ": done" << endl << endl;
		
	delete e1;
	delete tt1;
	delete tf1;
}

void chi2()
{
	string outchi = "mb_nce/results/chi2.txt";
	ofstream file(outchi.c_str());
	
	for (int i = 0; i < nof_ma; i++)
	{			
		stringstream temp;
		string x;
		
		temp << ma[i];
		temp >> x;
		
		string in = "mb_nce/results/rec_" + x + ".txt";
			
		double chi = calc_chi(in);

		file << ma[i] << " " << chi << endl;
		cout << ma[i] << " " << chi << endl;
		
		file << endl;
	}
	
	file.close();
}

double calc_chi(string in)
{
	double M[51][51] = {{0}};
	double E[51] = {0};
		
	read_rtab("mb_nce/tables/Mrev.txt", M);
	read_Erec(in, E);
		
	double res = 0;
	
	double dif[51];
		
	for (int i = 0; i < 51; i++)
		dif[i] = ::data[i] - E[i];
				
	for (int i = 0; i < 51; i++)
		for (int j = 0; j < 51; j++)
			res += dif[i]*M[i][j]*dif[j];
	
	return res;
}

void q2_calc(string rootC)
{
	string rootH = "mb_nce/root/H_qel.root";
	string txtH = rootH + ".txt";

	string txtC = rootC + ".txt";
	
	double H[50] = {0};
	double C[50] = {0};
	double CH2[50] = {0};
		
	double xsecH = crosssection(txtH);
	double xsecC = crosssection(txtC);
	
	if (waga)
	{
		xsecH = 1.0;
		xsecC = 1.0;
	}
	
	q2_calc(rootH, H);
	q2_calc(rootC, C);
	
	for (int i = 0; i < 50; i++)
	{
		H[i] *= 2.0 * xsecH / events / 0.04;
		C[i] *= 12.0 * xsecC / events / 0.04;
		CH2[i] = (H[i] + C[i]) / 14.0;
	}
	string q2out = "mb_nce/results/q2_dist.txt";
	ofstream file(q2out.c_str());
	
	for (int i = 0; i < 50; i++)
		file << i*0.04 << " " << CH2[i]*1e40 << endl <<
				(i+1)*0.04 << " " << CH2[i]*1e40 << endl;
	
	file.close();
}

void q2_calc(string in, double *res)
{
	TFile *tf1 = new TFile(in.c_str());
	TTree *tt1 = (TTree*)tf1->Get("treeout");
	event *e1   = new event();
		
	tt1->SetBranchAddress("e",&e1);
	
	for (int i = 0; i < events; i++)
	{
		tt1->GetEntry(i);
		
		if (e1->flag.qel == 1)
		{			
			double Q2 = -e1->q2()/1000000.0;
			int w = Q2/0.04;
			
			if (w < 50 and w >= 0)
			{
				if (waga)
					res[w] += e1->weight;
				else
					res[w]++;
			}
		}
		
		cout << in << ": " << 100*i/events << "%\r" << flush;
	}
	
	cout << in << ": done" << endl << endl;
		
	delete e1;
	delete tt1;
	delete tf1;
}

void sim_cc()
{
	string C = "mb_nce/root/C_cc.root";
	string H = "mb_nce/root/H_cc.root";
	
	if (noFile(C))
	{
		string command = "./bin/nuwro -o '" + C + "' " + parCC + fg + ccma + C;
		run(command);
	}
	
	if (noFile(H))
	{
		string command = "./bin/nuwro -o '" + H + "' " + parCC + fg + ccmaH + H;
		run(command);
	}
}

void ratio_calc(string rootCnc)
{
	string rootHnc = "mb_nce/root/H_qel.root";
	string txtHnc = rootHnc + ".txt";

	string txtCnc = rootCnc + ".txt";
	
	string rootHcc = "mb_nce/root/H_cc.root";
	string rootCcc = "mb_nce/root/C_cc.root";
	
	string txtHcc = rootHcc + ".txt";
	string txtCcc = rootCcc + ".txt";
	
	double Hnc[50] = {0};
	double Cnc[50] = {0};
	double Cnclike[50] = {0};
	double CH2nc[50] = {0};
	double CH2nclike[50] = {0};

	double Hcc[50] = {0};
	double Ccc[50] = {0};
	double Ccclike[50] = {0};
	double CH2cc[50] = {0};
	double CH2cclike[50] = {0};
	
	ratio_calcH(rootHnc, Hnc);
	ratio_calcH(rootHcc, Hcc);
	
	ratio_calcC(rootCnc, Cnc, Cnclike);
	ratio_calcC(rootCcc, Ccc, Ccclike);
	
	double xsecHnc = crosssection(txtHnc);
	double xsecCnc = crosssection(txtCnc);

	double xsecHcc = crosssection(txtHcc);
	double xsecCcc = crosssection(txtCcc);
	
	if (waga)
	{
		xsecHnc = 1.0;
		xsecCnc = 1.0;
		xsecHcc = 1.0;
		xsecCcc = 1.0;
	}
	
	for (int i = 0; i < 50; i++)
	{
		Hnc[i] *= 2 * xsecHnc / events / 0.04;
		Cnc[i] *= 12 * xsecCnc / events / 0.04;
		Cnclike[i] *= 12 * xsecCnc / events / 0.04;
		CH2nc[i] = (Hnc[i] + Cnc[i]) / 14;
		CH2nclike[i] = (Hnc[i] + Cnclike[i]) / 14;

		Hcc[i] *= 2 * xsecHcc / events / 0.04;
		Ccc[i] *= 12 * xsecCcc / events / 0.04;
		Ccclike[i] *= 12 * xsecCcc / events / 0.04;
		CH2cc[i] = (Hcc[i] + Ccc[i]) / 14;
		CH2cclike[i] = (Hcc[i] + Ccclike[i]) / 14;	
	}
	
	string rout = "mb_nce/results/ratio.txt";
	string routlike = "mb_nce/results/ratiolike.txt";
	
	ofstream file(rout.c_str());
	ofstream filelike(routlike.c_str());
	
	for (int i = 0; i < 50; i++)
	{
		file << i*0.04 << " " << (14.0/6.0)*CH2nc[i]/CH2cc[i] << endl << (i+1)*0.04 << " " << (14.0/6.0)*CH2nc[i]/CH2cc[i] << endl;
		filelike << i*0.04 << " " << (14.0/6.0)*CH2nclike[i]/CH2cclike[i] << endl << (i+1)*0.04 << " " << (14.0/6.0)*CH2nclike[i]/CH2cclike[i] << endl;
	}
	file.close();
	filelike.close();
}

void ratio_calcC(string in, double *res, double *reslike)
{
	TFile *tf1 = new TFile(in.c_str());
	TTree *tt1 = (TTree*)tf1->Get("treeout");
	event *e1   = new event();
		
	tt1->SetBranchAddress("e",&e1);
	
	for (int i = 0; i < events; i++)
	{
		tt1->GetEntry(i);
		
		bool nopi = (e1->nof(111) + e1->nof(211) + e1->nof(-211) == 0);
		bool nopilike = (e1->fof(111) + e1->fof(211) + e1->fof(-211) == 0);
		
		double Q2 = -e1->q2()/1000000.0;
		int w = Q2/0.04;
			
		if (w >= 50 or w < 0)
			continue;
		
		if (nopi)
		{
			if (waga)
				res[w] += e1->weight;
			else
				res[w]++;
		}
		
		if (nopilike)
		{
			if (waga)
				reslike[w] += e1->weight;
			else
				reslike[w]++;
		}
		
		cout << in << ": " << 100*i/events << "%\r" << flush;
	}
	
	cout << in << ": done" << endl << endl;
		
	delete e1;
	delete tt1;
	delete tf1;
}

void ratio_calcH(string in, double *res)
{
	TFile *tf1 = new TFile(in.c_str());
	TTree *tt1 = (TTree*)tf1->Get("treeout");
	event *e1   = new event();
		
	tt1->SetBranchAddress("e",&e1);
	
	for (int i = 0; i < events; i++)
	{
		tt1->GetEntry(i);
		
		bool nopi = (e1->nof(111) + e1->nof(211) + e1->nof(-211) == 0);
		
		double Q2 = -e1->q2()/1000000.0;
		int w = Q2/0.04;
			
		if (w >= 50 or w < 0)
			continue;
		
		if (nopi)
		{
			if (waga)
				res[w] += e1->weight;			
			else
				res[w]++;
		}
		
		cout << in << ": " << 100*i/events << "%\r" << flush;
	}
	
	cout << in << ": done" << endl << endl;
		
	delete e1;
	delete tt1;
	delete tf1;
}

void tm_sim(double ma)
{
	stringstream temp;
	string sma;
		
	temp << ma;
	temp >> sma;
			
	string HF = "mb_nce/root/H_qel.root";
	string resF = "mb_nce/root/C_res.root";
	string qelF = "mb_nce/root/C_qel_" + sma + ".root";
	
	if (noFile(HF))
	{
		string command = "./bin/nuwro -o '" + HF + "' " + parQEL + fg + ccma + ncma(1030) + H;
		run(command);
	}		
	
	if (noFile(resF))
	{
		string command = "./bin/nuwro -o '" + resF + "' " + parRES + fg + ccma + ncma(1350) + C + wf(7);
		run(command);
	}		
	
	if (noFile(qelF))
	{	
		string command = "./bin/nuwro -o '" + qelF + "' " + parQEL + fg + ccma + ncma(ma) + C + wf(7);
		run(command);
	}
}

void tm_calc(double ma)
{	
	stringstream temp;
	string x;
		
	temp << ma;
	temp >> x;

	string qelH = "mb_nce/root/H_qel.root";
	string qelHtxt = qelH + ".txt";

	string rdC = "mb_nce/root/C_res.root";
	string rdCtxt = rdC + ".txt";
	
	string qelC = "mb_nce/root/C_qel_" + x + ".root";
	string qelCtxt = qelC + ".txt";
	
	string out = "mb_nce/results/true_" + x + ".txt";
	string out2 = "mb_nce/results/rec_" + x + ".txt";
	
	if (noFile(out) or noFile(out2))
	{	
		double xsecH = crosssection(qelHtxt);
		double resxsec = crosssection(rdCtxt);
		double xsecC = crosssection(qelCtxt);
		
		if (waga)
		{
			xsecH = 1.0;
			xsecC = 1.0;
			resxsec = 1.0;
		}

		double result[5][51] = {{0}};

		calcH(qelH, result[0]);
		calcC(rdC, result);
		calcC(qelC, result);
			
		for (int i = 0; i < 51; i++)
		{
				result[4][i] *= resxsec/events/18.0;
				result[0][i] *= xsecH/events/18.0;
				result[1][i] *= xsecC/events/18.0;
				result[2][i] *= xsecC/events/18.0;
				result[3][i] *= xsecC/events/18.0;
		}
			
		norm(result[0], 2.0/14);
		norm(result[4], 12.0/14);
		norm(result[1], 12.0/14);
		norm(result[2], 12.0/14);
		norm(result[3], 12.0/14);
				
		double result2[51];
			
		true2rec(result, result2);
				
		ofstream file(out.c_str());
			
		for (int i = 0; i < 51; i++)
		{
			file << i*18 << " " << result[0][i] << " " << result[1][i] << " " << result[2][i] << " " << result[3][i] << " " << result[4][i] << endl;
			file << (i+1)*18 << " " << result[0][i] << " " << result[1][i] << " " << result[2][i] << " " << result[3][i] << " " << result[4][i] << endl;
		}
				
		file.close();

		ofstream file2(out2.c_str());
				
		for (int i = 0; i < 51; i++)
		{
			file2 << 40 + i*12 << " " << result2[i] << endl;
			file2 << 40 + (i+1)*12 << " " << result2[i] << endl;
		}

		file2.close();
	}
}

double tm_chi(double ma)
{
	stringstream temp;
	string x;
		
	temp << ma;
	temp >> x;
		
	string in = "mb_nce/results/rec_" + x + ".txt";

	double chi = calc_chi(in);
								
	return chi;
}

void norm(double *tab, double x)
{	
	double factor = deltaT * Nn * POT * flux;

	for (int j = 0; j < 51; j++)
			tab[j] *= factor*x;
}

void normHE(double *tab, double x)
{	
	double factor = deltaT * Nn * POT * flux;

	for (int j = 0; j < 30; j++)
			tab[j] *= factor*x;
}

void run(string com) //run external program "com"
{
	FILE *fp;
	fp = popen(com.c_str(), "w");
	pclose(fp);
}

double crosssection (string filename) //read xsec from *.root.txt
{
	string y;
	double x;
		
	vector<double> v(8);
	ifstream Input (filename.c_str());

	if (Input)
	{
		do
		{
			getline (Input, y);
			for (int j = 0; j < 8; j++)
				{
					for (int k = 0; k < 3; k++)
					{
						Input>>y;
					}
					Input>>v[j];
				}
		}while (Input);
	}
	
	double result = 0;
	
	for (int i = 0; i < 8; i++)
	{
		result += v[i];
	}
	
	return result;
}

string ncma(double val) //create input string from double for Ma
{
	stringstream temp;
	string x;
	
	temp << val;
	temp >> x;
	
	string ncma = "-p 'qel_nc_axial_mass = " + x + "' -p 'qel_s_axial_mass = " + x + "' ";
	return ncma;
}

string wf(double val) //create input string from double for work function
{
	stringstream temp;
	string x;
	
	temp << val;
	temp >> x;
	
	string ncma = "-p 'kaskada_w = " + x + "' ";
	return ncma;
}

string ds(double val) //create input string from double for delta_s
{
	stringstream temp;
	string x;
	
	temp << val;
	temp >> x;
	
	string ncma = "-p 'delta_s = " + x + "' ";
	return ncma;
}

bool noFile(string filename) //returns true if a file does not exist
{
	fstream plik;
	plik.open(filename.c_str());
	
	if (plik.is_open())
	{
		plik.close();
		return false;
	}
	
	plik.close();
	return true;	
}
