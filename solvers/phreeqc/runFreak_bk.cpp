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
	char sep ='/';
	std::ifstream inputInit{cur_dir+sep+"phqinit.txt"};
	std::vector<int> ph_init{std::istream_iterator<int>{inputInit}, {}}; 
	//nxyz=ph_init[0];ph_ncomp=ph_init[1];ph_gcomp=ph_init[2];ph_nsolu=ph_init[3];
	nxyz=2;ph_ncomp=8;ph_gcomp=0;ph_nsolu=2;
	std::vector<int> data={2,8,0,2,1, -1,0,1, 0,-1,0,-1,0,-1,0,-1,0,-1,0,1};

	/*
	freak.mix=true;
	std::ifstream inputMixs{cur_dir+sep+"phqmixs.txt"};
	std::vector<int> ph_mixs{std::istream_iterator<int>{inputMixs}, {}}; 
	freak.setMixs(ph_mixs);
	std::ifstream inputMixf{cur_dir+sep+"phqmixf.txt"};
	std::vector<double> ph_mixf{std::istream_iterator<double>{inputMixf}, {}}; 
	freak.setMixf(ph_mixf);*/
	freak.setDB(cur_dir+sep+"phreeqc.dat");
	//make a first init to calculate the inital solutions (poro required for the initial solutions)
	freak.setChemFile(cur_dir+sep+"initChem.pqi"); //Info << "initCh read " << endl;
	freak.setData(data);
	std::vector<double> poro(ph_nsolu,0.25);
	std::vector<double> wsat(ph_nsolu,0.9999);
	freak.setPoro(poro);
	freak.setWsat(wsat);
	freak.init(); 
	
	//first step make 50 time steps for solution 1
	int icell,it,i,j,k;int nsp=ph_ncomp-4;
	int nt=50;std::vector<double> c(nt*nsp,0.);
	std::vector<double> cmax(nsp,0.);
	for (it=1;it<nt;it++) { //50 time steps
		freak.setTstep(it*86400.);
		freak.run();icell=1;
		for (i=0;i<nsp;i++) { 
			c[it*nsp+i]=freak.c[(i+4)*ph_nsolu+icell];
			cmax[i] = std::max(cmax[i],c[it*nsp+i]);
			}
	}
	std::vector<double> sbase0(4);
	for (i=0;i<4;i++) {sbase0[i] = freak.c[i*nxyz+0];}
	// then from these 50 solutions, mix them with solution zero by 1/10
	std::vector<double> c0(nsp,0.);int nk=11;
	for (i=0;i<nsp;i++) {c0[i]=freak.c[(i+4)*nxyz+0];} // intial solution
	std::vector<double> c1(nt*nsp*nk,0.);
	for (it=1;it<nt;it++) {
		for (k=0;k<nk;k++) {
			for (i=0;i<nsp;i++) {
				c1[it*nk*nsp+k*nsp+i]=k*0.1*c0[i]+(10-k)*0.1*c[it*nsp+i];
				}
		}
	}

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