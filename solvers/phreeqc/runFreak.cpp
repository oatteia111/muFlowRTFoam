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
	std::vector<double> poro(ph_nsolu,0.25);std::vector<double> t_ph(ph_nsolu,20.);
	std::vector<double> wsat(ph_nsolu,0.9999);
	freak.setPoro(poro);freak.setTemp(t_ph);
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
	int ntest=10; nxyz = nmix*ntest;double x,r; //for each mix we will have 5 different options of use of phase, exch, surf
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
	
	//std::ifstream inputMNmnx{cur_dir+sep+"nminx.txt"};
	//std::vector<int> nminx{std::istream_iterator<int>{inputNminx}, {}}; //nb or equilibirm phases in selected output
	std::ofstream nnData(cur_dir+sep+"nnData.txt");
	std::ofstream nnTarget(cur_dir+sep+"nnTarget.txt");
	//run phreeqc
	std::vector<double> temp(nxyz,25.);
	//freak.setTemp(temp);
	//first add concentrations and make a first run to equilibrate
	std::vector<double>c_ph(nxyz*ph_ncomp);
	for (j=0;j<nmix;j++)
		{
		for (k=0;k<ntest;k++)
			{
			j1 = j*ntest+k;r=0.9+std::rand()/4e10;
			for (i=0;i<4;i++) {c_ph[i*nxyz+j1] = c[j*ph_ncomp+i];}
			for (i=4;i<ph_ncomp;i++) {c_ph[i*nxyz+j1] = c[j*ph_ncomp+i]*r;}
			}
		} std::cout<<"c_ph filled \n";
	freak.setC(c_ph);
	freak.run();std::cout<<"end run 1\n";
	//reading selected output
	freak.getSelOutput();int nse1=freak.nselect-1;
	std::cout<<" nse1 "<<nse1<<" "; //base
	/*std::vector<double>s_ph1(nxyz*nse1);
	for (j=0;j<nxyz;j++)
		{
		s_ph1[j]=freak.spc[j];//pH, we remove pe
		for (i=1;i<nse1;i++) {s_ph1[i*nxyz+j]=freak.spc[(i+1)*nxyz+j];}
		}
	// slightly modify the concentrations, equilibrate, to see the variations of c and solids, data contains c
	for (j=0;j<nxyz;j++)
		{
		r=0.99+std::rand()/2e11;
		for (i=4;i<ph_ncomp;i++) {x = std::max(c_ph[i*nxyz+j]*r,1e-15);c_ph[i*nxyz+j] = x;nnData<<x<<" "; }
		for (i=0;i<nse1;i++) {nnData<<s_ph1[i*nxyz+j]<<" ";} //selected out species
		nnData<<"\n";
		}
	std::cout<<"end sel out \n";*/
	//freak.setPoro(poro);
	//freak.setWsat(wsat);
	freak.setC(c_ph);
	//freak.setTemp(temp);
	freak.run();
	std::cout<<"end run 2 \n";
	/*freak.getSelOutput();std::vector<double>s_dff1(nxyz*nse1); //base
	for (j=0;j<nxyz;j++)
		{
		s_dff1[j]=freak.spc[j]-s_ph1[j]; //pH, we remove pe
		for (i=1;i<nse1;i++) {s_dff1[i*nxyz+j]=freak.spc[(i+1)*nxyz+j]-s_ph1[i*nxyz+j];}
		}
	//now write the conc differences in target
	for (j=0;j<nxyz;j++)
		{
		for (i=4;i<ph_ncomp;i++) {nnTarget<<freak.c[i*nxyz+j]-c_ph[i*nxyz+j]<<" ";}//ntest times the same solution
		for (i=0;i<nse1;i++) {nnTarget<<s_dff1[i*nxyz+j]<<" ";}
		//for (i=0;i<nse2;i++) {nnTarget<<s_dff2[i*nxyz+j]<<" ";}
		nnTarget<<"\n";
		}
	std::cout<<"stored in target \n";*/
}
