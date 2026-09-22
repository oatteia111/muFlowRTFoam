#if defined(USE_MPI)
#include <mpi.h>
#endif

#ifndef initPhreeqc_H
#define initPhreeqc_H

#include <stdlib.h>
#include <vector>
#include <fstream>
#include <iterator>
#include <sstream>
#include <string>
#include <chrono> 

// no includes here for windows bug

class PhreeqcRM; // Déclaration anticipée magique pour Windows

//int nthreads = 3;
int nxyz,ph_ncomp,ph_nspc,ph_gcomp,ph_nsolu;

////////////////////// start phreeqc
class my_phq
{
public :
    my_phq(); // Juste un point-virgule, pas d'accolades avec du code

	PhreeqcRM *PhreeqcRM_ptr;
#ifdef USE_MPI
	MPI_Comm rm_commxx;
#endif

	std::string DB;
	std::vector<int> data;
	std::string ChemFile;
	std::vector<double> c;
	std::vector<double> bc_conc;
	std::vector<double> bc_species;
	std::vector<double> g;
	std::vector<double> gm;
	std::vector<double> temp;
	std::vector<double> p;
	std::vector<double> spc;
	double tstep=1;
	std::vector<double> poro;
	std::vector<double> wsat;
	std::vector<double> gvol;
	std::vector<double> foc;
	std::vector<std::string> comp;
	std::vector<std::string> gcomp;
	int iGwater=-1;
	bool EK;
	int nspc,nselect,nphase,nexch,nsurf;
	std::vector<double> diff25,deltaMvol;
	std::vector<double> z;
	std::vector<std::string> spcNames;

	void setDB(const std::string& db);
	void setData(std::vector<int> data){this->data = data;}
	void setChemFile(std::string chfile){this->ChemFile = chfile;}
	void setC(std::vector<double> c){this->c = c;}
	void setG(std::vector<double> g){this->g = g;}
	void setGm(std::vector<double> gm){this->gm = gm;}
	void setTemp(std::vector<double> temp){this->temp = temp;}
	void setP(std::vector<double> p){this->p = p;}
	void setSpc(std::vector<double> spc){this->spc = spc;}
	void setTstep(double tstep){this->tstep = tstep;}
	void setPoro(std::vector<double> poro){this->poro = poro;}
	void setWsat(std::vector<double> wsat){this->wsat = wsat;}
	void setGvol(std::vector<double> gvol){this->gvol = gvol;}
	void setFoc(std::vector<double> foc){this->foc = foc;}
	void setNwat(int){this->iGwater = iGwater;}

    int phqInit(my_phq & freak);
	int phqRun(my_phq & freak);
	int getSelOutput(my_phq & freak);
	int phqTest(std::string time) ;
	int phqRestart(my_phq freak, std::vector<int> & ph_data, std::string time,std::vector<double> c_ph,std::vector<double> g_ph) ;
	void register_basic_callback(void *cookie);
	static double my_basic_callback(double x1, double x2, const char *str, void *cookie);

private:
    PhreeqcRM* phreeqc_instance; // Si vous utilisez un pointeur, c'est direct.
};
#endif