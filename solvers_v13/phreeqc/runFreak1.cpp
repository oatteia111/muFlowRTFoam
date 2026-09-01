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
int i,j;
#include "initPhreeqc.H"
my_phq freak;

int main(int argc, char *argv[])
{
	char sep ='/';std::cout<<"start \n";
	std::ifstream inputInit{cur_dir+sep+"phqinit.txt"};
	std::vector<int> ph_init{std::istream_iterator<int>{inputInit}, {}}; 
	nxyz=10;ph_nsolu=10;
	//nxyz=ph_init[0]=nxyz;
	ph_ncomp=ph_init[1];ph_gcomp=ph_init[2];//ph_nsolu=ph_init[3];
	std::vector<int> data={nxyz,ph_ncomp,ph_gcomp,ph_nsolu,1, 0,0, 0,-1, 0,-1, 0,-1, 0,0, 0,-1, 0,-1};

	freak.setDB(cur_dir+sep+"phreeqc.dat");
	//make a first init to calculate the inital solutions (poro required for the initial solutions)
	freak.setChemFile(cur_dir+sep+"initChem.pqi"); //Info << "initCh read " << endl;
	freak.setData(data);
	std::vector<double> poro(ph_nsolu,0.25);std::vector<double> t_ph(ph_nsolu,0.);
	for (i=0;i<nxyz;i++) {t_ph[i]=10.5+i*2.;}
	std::vector<double> wsat(ph_nsolu,.9999);
	freak.setPoro(poro);freak.setTemp(t_ph);
	freak.setWsat(wsat);
	freak.init(); 
	std::vector<double> gvol(ph_nsolu,0.24);
	std::vector<double> gm_ph(ph_nsolu*ph_gcomp,0.01);
	freak.setGvol(gvol);freak.setGm(gm_ph);freak.setTstep(100.);freak.setTemp(t_ph);
	freak.run();
	for (j=0;j<nxyz;j++) { for (i=0;i<2;i++) {std::cout<<"g "<<freak.g[i*nxyz+j]<<" gm "<<freak.gm[i*nxyz+j]<<" ";} std::cout<<"\n";}
	std::cout<<"end run \n";
}
