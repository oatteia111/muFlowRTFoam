#include <stdlib.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <iterator>
#include <sstream>
#include <string>

#include <algorithm>
#include <iterator>
#include <random>

//////////////////// find local dir
#include <unistd.h>
#define GetCurrentDir getcwd

std::string get_current_dir() {
   char buff[FILENAME_MAX]; //create string buffer to hold path
   GetCurrentDir( buff, FILENAME_MAX );
   std::string current_working_dir(buff);
   return current_working_dir;
}
std::string cur_dir = get_current_dir();

#include "initPhreeqc.H"
my_phq freak;

int main(int argc, char *argv[])
{
	char sep ='/';std::cout<<"start \n";
	std::ifstream inputInit{cur_dir+sep+"phqinit.txt"};
	std::vector<int> ph_init{std::istream_iterator<int>{inputInit}, {}}; 
	nxyz=ph_init[0];ph_ncomp=ph_init[1];ph_gcomp=ph_init[2];ph_nsolu=ph_init[3];
	//std::vector<int> data={nxyz,ph_ncomp,0,ph_nsolu,1, -1,0,1, 0,-1,0,-1,0,-1,0,-1,0,-1,0,1};

	freak.setDB(cur_dir+sep+"phreeqc.dat");
	//make a first init to calculate the inital solutions (poro required for the initial solutions)
	freak.setChemFile(cur_dir+sep+"initChem.pqi"); //Info << "initCh read " << endl;
	freak.setData(ph_init);
	std::vector<double> poro(ph_nsolu,0.25);
	std::vector<double> wsat(ph_nsolu,0.9999);
	freak.setPoro(poro);
	freak.setWsat(wsat);
	freak.init(); 
	int i,j,k;int nmix=100;
	std::vector<double> c(nmix*ph_ncomp,0.);
	for (k=0;k<ph_nsolu;k++) //k is the solution number, we keep the orignal solutions with no mix
		{	
		//std::cout<<" k "<<k<<" ";
		for (i=0;i<ph_ncomp;i++) {c[k*ph_ncomp+i]= freak.c[i*ph_nsolu+k];} //std::cout<<c[k*ph_ncomp+i]<<" ";
		//std::cout<<"\n";
		}
	
	//create random 100 solutions mixed from the n initial solutions
	// the proportions should be around 1/nsolu where nsolu is the nb of solutions
	int icell,it;//int nsp=ph_ncomp-4;
	float prop,totprop;
	for (j=ph_nsolu;j<nmix;j++) // the first solution are already stored
		{
		totprop=0;
		for (k=0;k<ph_nsolu;k++) //k is the solution number
			{	
			prop = std::rand()/2e9/(ph_nsolu-1); totprop += prop; //std::cout<<" j "<<j <<" k "<<k<<" prop "<<prop<<"\n";
			if (totprop>1.) {prop -= totprop-1;}
			if ((k==ph_nsolu-1)&&(totprop<=1.)) {prop += 1-totprop;}
			for (i=0;i<ph_ncomp;i++) {c[j*ph_ncomp+i] += prop*freak.c[i*ph_nsolu+k];} // in c one line=one solu (!=phreeqc)
			}
		/*std::cout<<" j "<<j<<" ";
		for (i=0;i<ph_ncomp;i++) {std::cout<<c[j*ph_ncomp+i]<<" ";}
		std::cout<<"\n";*/
		}
	std::cout<<"props done \n";
	// now create input with mixe sin contact to combination of phase, exchange and surf
	int nphase=6;int nsurf=6;int nexch=6;int j1,ip,ie,is; // nb in example trip these phases... are all the same, but can be different
	int ntest=5; nxyz = nmix*ntest;double x,r; //for each mix we will have 5 different options of use of phase, exch, surf
	std::vector<int> data={nxyz,ph_ncomp,0,nxyz,1, 0, 0}; //up to now we set all to soutions to 0
	std::vector<int> phase(nxyz),exch(nxyz),surf(nxyz);
	// 
	for (j=0;j<nmix;j++)
		{
		for (k=0;k<ntest;k++)
			{
			ip = std::rand()%nphase; phase.push_back(ip);
			ie = std::rand()%nexch; exch.push_back(ie);
			is = std::rand()%nsurf; surf.push_back(is);
			}
		}
	data.push_back(-1); //phases
	for (j=0;j<nxyz;j++) {data.push_back(phase[j]);}
	data.push_back(-1); //exchanger
	for (j=0;j<nxyz;j++) {data.push_back(exch[j]);}
	data.push_back(-1); //surfaces
	for (j=0;j<nxyz;j++) {data.push_back(surf[j]);}
	std::vector<int> dend{0,-1,0,-1,0,-1}; // gases solid slution, kinetics
	data.insert(data.end(), dend.begin(), dend.end());std::cout<<"solids done \n";
	freak.setData(data);
	poro.resize(nxyz,0.25);
	wsat.resize(nxyz,0.9999);
	freak.setPoro(poro);
	freak.setWsat(wsat);
	freak.init(); std::cout << "freak H2O_0 " << freak.c[0] << "\n";
	
	std::ofstream nnData(cur_dir+sep+"nnData.txt");
	std::ofstream nnTarget(cur_dir+sep+"nnTarget.txt");
	//run phreeqc
	std::vector<double> temp(nxyz,25.);
	freak.setTemp(temp);
	//first add concentrations and make a first run to equilibrate
	std::vector<double>c_ph(nxyz*ph_ncomp);
	for (j=0;j<nmix;j++)
		{
		for (k=0;k<ntest;k++)
			{
			j1 = j*ntest+k;r=0.99+std::rand()/2e11;
			for (i=0;i<ph_ncomp;i++) {c_ph[i*nxyz+j1] = c[j*ph_ncomp+i]*r;}//ntest times the same solution
			}
		} std::cout<<"c_ph filled \n";
	freak.setC(c_ph);
	freak.run();std::cout<<"end run 1\n";
	// slightly modify the concentrations, equilibrate, to see the variations of solids, data contains the delta in c
	std::vector<double>s_ph(nxyz*freak.nselect);
	for (j=0;j<nxyz;j++)
		{
		r=0.99+std::rand()/2e11;
		for (i=0;i<ph_ncomp;i++) {x = c_ph[i*nxyz+j];c_ph[i*nxyz+j] = x*r;nnData<<x*r<<" ";}
		for (i=0;i<freak.nselect;i++) {s_ph[i*nxyz+j]=freak.spc[i*nxyz+j];nnData<<s_ph[i*nxyz+j]<<" ";}
		nnData<<"\n";
		}
	freak.setC(c_ph);
	freak.run();std::cout<<"rn 2 \n";
	//now write the results in target (conc differences)
	for (j=0;j<nxyz;j++)
		{
		for (i=0;i<ph_ncomp;i++) {nnTarget<<freak.c[i*nxyz+j]-c_ph[i*nxyz+j]<<" ";}//ntest times the same solution
		for (i=0;i<freak.nselect;i++) {nnTarget<<freak.spc[i*nxyz+j]-s_ph[i*nxyz+j]<<" ";}
		nnTarget<<"\n";
		}
}
	
/*
	// use these concentrations to make a run for a given time step and get dC/dt
	// make conc relative, here we keep c0 as cmin, and store in nndata
	std::ofstream nnData(cur_dir+sep+"nnData.txt");
	std::ofstream nnTarget(cur_dir+sep+"nnTarget.txt");
	nxyz = nt*nk;
	data={nxyz,8,0,nxyz,1, 0,0, 0,-1,0,-1,0,-1,0,-1,0,-1, 0,1};
	freak.setData(data);
	poro.resize(nxyz,0.25);
	wsat.resize(nxyz,0.9999);
	freak.setPoro(poro);
	freak.setWsat(wsat);
	freak.init(); std::cout << "freak H2O_0 " << freak.c[0] << "\n";
	double dt = 86400.;
	freak.setTstep(dt);
	std::vector<double> c_ph(nxyz*ph_ncomp,0);
	for (int j=0;j<nxyz;j++) {
		for (i=0;i<4;i++) {c_ph[i*nxyz+j] = freak.c[i*nxyz+j];}
		for (i=0;i<nsp;i++) {c_ph[(i+4)*nxyz+j] = c1[j*nsp+i];}
		//if (i==0) {std::cout<<c_ph[0]<<"\n";}
		}
	freak.setC(c_ph);
	freak.run();std::cout<<"end run \n";
//shuffle
std::vector<int> idx(nt*nk);
std::iota(idx.begin(), idx.end(), 0);		
std::random_device rd;
std::mt19937 g(rd());
std::shuffle(idx.begin(), idx.end(), g);
	std::vector<double> dC(nxyz*nsp,0.);
	for (int j=1;j<nxyz;j++) {
		int j1 = idx[j];
		for (i=0;i<nsp;i++) {
			nnData << (c1[j1*nsp+i]-c0[i])/(cmax[i]-c0[i]) << " "; //crelative
			dC[j1*nsp+i]=(freak.c[(i+4)*nxyz+j1]-c1[j1*nsp+i]);
			nnTarget << dC[j1*nsp+i]/dt*1e10 << " ";
		}
		nnData << "\n";nnTarget << "\n";
	}
}
*/