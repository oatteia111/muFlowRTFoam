
#if defined(USE_MPI)
#include <mpi.h>
#endif
#include <stdlib.h>
#include <iostream>
#include <string>
#include <vector>

#include <fstream>
#include <iterator>
#include <sstream>

#include "PhreeqcRM.h"
#include "IPhreeqc.hpp"
#include "IPhreeqcPhast.h"
int nthreads = 3;

int nxyz,ph_ncomp,ph_gcomp,ph_nsolu,sunits;
double my_basic_callback(double x1, double x2, const char *str, void *cookie);
void register_basic_callback(void *cookie);
										  																		  										 
class my_phq
{
public:
	PhreeqcRM* PhreeqcRM_ptr;
#ifdef USE_MPI
	MPI_Comm rm_commxx;
#endif
	std::string DB;
	void setDB(std::string db) { this->DB = db; }
	std::vector<int> data;
	void setData(std::vector<int> data) { this->data = data; }
	std::string ChemFile;
	void setChemFile(std::string chfile) { this->ChemFile = chfile; }
	std::vector<double> c0;
	void setC0(std::vector<double> c0) { this->c0 = c0; }
	std::vector<double> c;
	void setC(std::vector<double> c) { this->c = c; }
	std::vector<double> g;
	void setG(std::vector<double> g) { this->g = g; }
	std::vector<double> gm;
	void setGm(std::vector<double> gm) { this->gm = gm; }
	std::vector<double> temp;
	void setTemp(std::vector<double> temp) { this->temp = temp; }
	std::vector<double> p;
	void setP(std::vector<double> p) { this->p = p; }
	std::vector<double> spc;
	void setSpc(std::vector<double> spc) { this->spc = spc; }
	//std::vector<double> solu_conc;
	//void setSconc(std::vector<double> solu_conc){this->solu_conc = solu_conc;}
	double tstep = 1;
	void setTstep(double tstep) { this->tstep = tstep; }
	std::vector<double> rv;
	void setRV(std::vector<double> rv) { this->rv = rv; }
	std::vector<double> poro;
	void setPoro(std::vector<double> poro) { this->poro = poro; }
	std::vector<double> wsat;
	void setWsat(std::vector<double> wsat) { this->wsat = wsat; }
	std::vector<double> gvol;
	void setGvol(std::vector<double> gvol) { this->gvol = gvol; }
	std::vector<std::string> comp;
	std::vector<std::string> gcomp;						   							
	int iGwater = -1;
	void setNwat(int) { this->iGwater = iGwater; }
	bool EK;
	int nspc, nselect;
	std::vector<double> diff25;
	std::vector<double> z;
	std::vector<std::string> spcNames;
};

int phqInit(my_phq & freak)
{
	// !!! the first nxyz cells are true cells while the ph_nsolu following ones are for solutions
	int nxyz = 3;int nthreads = 2;
	int i, j, ph_ncomp, ph_gcomp, ph_nsolu;
	std::vector<int> ph_init = freak.data;
	nxyz = ph_init[0];ph_ncomp = ph_init[1];ph_gcomp = ph_init[2];ph_nsolu = ph_init[3];
	int sunits = ph_init[4];
	PhreeqcRM* phreeqc_rm = freak.PhreeqcRM_ptr;
	std::cout << "start of phqInit " << nxyz << std::endl;
	try
	{
		IRM_RESULT status;
		// Set properties
		status = phreeqc_rm->SetErrorHandlerMode(1);
		status = phreeqc_rm->SetComponentH2O(true);
		status = phreeqc_rm->SetRebalanceFraction(0.0);
		status = phreeqc_rm->SetRebalanceByCell(true);
		phreeqc_rm->UseSolutionDensityVolume(false);
		phreeqc_rm->SetPartitionUZSolids(false);
		std::cerr << "before  SetFilePrefix " << std::endl;
		// Open files
		status = phreeqc_rm->SetFilePrefix("phq1_cpp");
		//phreeqc_rm->OpenFiles();
		// Set concentration units
		status = phreeqc_rm->SetUnitsSolution(2);           // 1, mg/L; 2, mol/L; 3, kg/kgs
		status = phreeqc_rm->SetUnitsPPassemblage(sunits);       // 0, mol/L cell; 1, mol/L water; 2 mol/L rock
		status = phreeqc_rm->SetUnitsExchange(sunits);           // 0, mol/L cell; 1, mol/L water; 2 mol/L rock
		status = phreeqc_rm->SetUnitsSurface(sunits);            // 0, mol/L cell; 1, mol/L water; 2 mol/L rock
		status = phreeqc_rm->SetUnitsGasPhase(1);           // 0, mol/L cell; 1, mol/L water; 2 mol/L rock
		status = phreeqc_rm->SetUnitsSSassemblage(sunits);       // 0, mol/L cell; 1, mol/L water; 2 mol/L rock
		status = phreeqc_rm->SetUnitsKinetics(1);           // 0, mol/L cell; 1, mol/L water; 2 mol/L rock

		if (freak.EK) {status = phreeqc_rm->SetSpeciesSaveOn(true);}

		status = phreeqc_rm->SetRepresentativeVolume(freak.rv);   // RV is one dm3 of medium
		std::cerr << "before  SetPorosity " << freak.poro.size() << std::endl;
		status = phreeqc_rm->SetPorosity(freak.poro);
		std::cout << "end status \n";
		// Open for writing file
		status = phreeqc_rm->SetFilePrefix("phq1");
		//phreeqc_rm->OpenFiles();
		
		//load database and chemistry files and make the first run of the chem file
		status = phreeqc_rm->SetPrintChemistryOn(false, false, false); // workers, initial_phreeqc, utility
		status = phreeqc_rm->LoadDatabase(freak.DB);
		std::cout << "dbase opened " << "\n";
		bool workers = true;             // Worker instances do the reaction calculations for transport
		bool initial_phreeqc = true;     // InitialPhreeqc instance accumulates initial and boundary conditions
		bool utility = true;             // Utility instance is available for processing
		status = phreeqc_rm->RunFile(workers, initial_phreeqc, utility, freak.ChemFile);
		std::cerr << "After RunFile." << std::endl;

		int ph_ncomp = phreeqc_rm->FindComponents();freak.data[1]=ph_ncomp; // NECESSARY don't remove
		const std::vector<std::string>& gcompn = phreeqc_rm->GetGasComponents(); //std::cout<<"phq gcomp0 "<<gcomp[0]<<"\n";
		const std::vector<std::string>& compn = phreeqc_rm->GetComponents();
		std::cout << "phq comp "; for (i = 0;i < ph_ncomp;i++) { std::cout << compn[i] << " "; } std::cout << "\n";

		// Set array of initial conditions by reading the solution... numbers in data
		std::vector<int> ic1; //ic1 for the whole domain + iniitla solutions at the end
		ic1.resize(nxyz* 7, -1);std::cout<<"ic1 size "<<ic1.size()<<"\n";
		int icnt = 4; // to consider the first 5 numbers : nxyz,ncomp,gcomp,nsolu,units (++ in the first line of the for makes 5 numbers)
		
		for (i = 0; i < 7; i++)  //Solution,Equilibrium phases,Exchange,Surface,Gas phase,Solid solutions,Kinetics
		{
			icnt++;
			int test = freak.data[icnt];//std::cout<<"index "<<i<<" " <<test<<"\n";
			if (test == 0)  // if 0 just one value to read
			{icnt++;
			for (j=0;j<nxyz;j++) {ic1[i*nxyz+j] = freak.data[icnt];}
			}
			else // several values
			{
			for (j=0;j<nxyz;j++) {icnt++;ic1[i*nxyz+j] = freak.data[icnt];}
			}
		} 

		status = phreeqc_rm->InitialPhreeqc2Module(ic1);
		std::cout<<"end fill ic1 \n";
		int ncell = phreeqc_rm->GetChemistryCellCount();
		std::cout << "ncell " << ncell << " ph_ncomp " << ph_ncomp << " starting init phreeqc \n";
		//phreeqc_rm->SetPrintChemistryOn(true,false,false);
		//status = phreeqc_rm->SetTemperature(freak.temp);
		status = phreeqc_rm->RunCells();
		std::cout << " cells run \n";
		
		freak.wsat.resize(nxyz);
		status = phreeqc_rm->GetSaturation(freak.wsat);
		assert(compn.size() == (size_t)phreeqc_rm->GetComponentCount());
		std::vector<int> v(ph_nsolu);
		for (i=0;i<ph_nsolu;i++) {v[i]=i;}
		std::vector<double>c0(ph_nsolu * compn.size());
		status = phreeqc_rm->InitialPhreeqc2Concentrations(c0,v);
		freak.setC0(c0);
		
		std::vector<double>c(nxyz * compn.size());
		status = phreeqc_rm->GetConcentrations(c);
		freak.setC(c);

		if (ph_gcomp > 0) {
			status = phreeqc_rm->GetSaturation(freak.wsat);
			status = phreeqc_rm->GetGasCompPressures(freak.g);
			status = phreeqc_rm->GetGasCompMoles(freak.gm);
			status = phreeqc_rm->GetGasPhaseVolume(freak.gvol);
			std::cout << "in phq poro " << freak.poro[1] << "\n";
			for (size_t i = 0; i < compn.size(); i++)
			{
				std::cerr << compn[i] << " " << freak.c[i * nxyz + 1] << std::endl;
			}
			//std::cout<<"gvol 0 "<<freak.gvol[0]<<"gvol 1 "<<freak.gvol[1]<<"\n";
			std::cout<<"g 0 ";for (int i=0;i<ph_gcomp;i++) {std::cout <<freak.g[i*nxyz+0]<<" ";} std::cout<<"\n";
			std::cout<<"g 1 ";for (int i=0;i<ph_gcomp;i++) {std::cout <<freak.g[i*nxyz+1]<<" ";} std::cout<<"\n";
			for (size_t i = 0;i < gcompn.size();i++) {
				if (gcompn[i] == "H2O(g)") { freak.iGwater = i; }
			}
		}
		std::cout << "c ";for (int i=0; i<ph_ncomp;i++){std::cout<< freak.c[i*nxyz+0]<<" ";} std::cout<<"\n"; //ic1 seems ok
		std::cout << "end phq init " << "\n";

		// Get pointer to worker
		/*int iphreeqc_result;
		const std::vector<IPhreeqcPhast*> w = phreeqc_rm->GetWorkers();
		w[0]->AccumulateLine("Delete; -all");
		iphreeqc_result = w[0]->RunAccumulated();
		status = phreeqc_rm->CloseFiles();*/
	}
	catch (PhreeqcRMStop)
	{
		std::string e_string = "runf failed with an error in PhreeqcRM.";
		std::cerr << e_string << std::endl;
	#ifdef USE_MPI
		MPI_Abort(MPI_COMM_WORLD, 1);
	#endif
		return IRM_FAIL;
	}
	catch (...)
	{
		std::string e_string = "runf failed with an unhandled exception.";
		std::cerr << e_string << std::endl;
	#ifdef USE_MPI
		MPI_Abort(MPI_COMM_WORLD, 1);
	#endif
		return IRM_FAIL;
	}
	return EXIT_SUCCESS;
};

	//--------------------------- RUN ------------------------------
	
int phqRun(my_phq & freak)
	{
	try
	{
		IRM_RESULT status;
		int nxyz = freak.PhreeqcRM_ptr->GetGridCellCount();
		//int ph_ncomp = freak.PhreeqcRM_ptr->FindComponents();
		int ph_ncomp = freak.data[1];
		//freak.PhreeqcRM_ptr->OpenFiles();
		const std::vector<std::string>& compn = freak.PhreeqcRM_ptr->GetComponents();
		const std::vector<std::string>& gcompn = freak.PhreeqcRM_ptr->GetGasComponents();
		int ph_gcomp = (int)gcompn.size();
		std::cout << "phq ncomp sz " << ph_ncomp << " gcomp " << gcompn.size() << " ws " << freak.wsat[0] << "\n";
		std::cout << "phq comp     "; for (size_t i = 0;i < ph_ncomp;i++) { std::cout << compn[i] << " "; } std::cout << "\n";
		//int nxyz = freak.data[0];int ph_gcomp = freak.data[2];int ph_nsolu = freak.data[3];int nstot = nxyz + ph_ncomp;
		std::cout << "in phq nxyz " << nxyz << " size c " << freak.c.size() << " poro " << freak.poro[0] << std::endl;
		std::cerr << "Temp sz: " << freak.temp.size() << std::endl;
		
		status = freak.PhreeqcRM_ptr->SetTemperature(freak.temp);std::cout << "temp done ";
		if (gcompn.size() > 0) {
			std::cout << "nb gcomp before " << gcompn.size() << std::endl;//<<" gvol 3 "<<freak.gvol[0]<<" sat "<<freak.wsat[0]<<"\n";
			//std::cout<<"g 0 ";for (int i=0;i<ph_gcomp;i++) {std::cout <<freak.g[i*nxyz+0]<<" ";} std::cout<<"\n";
			//std::cout<<"gm 0 ";for (int i=0;i<ph_gcomp;i++) {std::cout <<freak.gm[i*nxyz+0]<<" ";} std::cout<<"\n";
			status = freak.PhreeqcRM_ptr->SetSaturation(freak.wsat); // keep this order to set values
			status = freak.PhreeqcRM_ptr->SetGasPhaseVolume(freak.gvol);
			status = freak.PhreeqcRM_ptr->SetGasCompMoles(freak.gm);
			//status = freak.PhreeqcRM_ptr->SetPressure(freak.p); // does not seem to be usefull 
		}
		std::cout << "EK " << freak.EK << "\n";
		if (freak.EK) 
		{ status = freak.PhreeqcRM_ptr->SpeciesConcentrations2Module(freak.spc); }
		else 
		{ status = freak.PhreeqcRM_ptr->SetConcentrations(freak.c); }
	
		status = freak.PhreeqcRM_ptr->SetTimeStep(freak.tstep);std::cout << " tstep " << freak.tstep << "\n";

		//----- run here ------
		status = freak.PhreeqcRM_ptr->RunCells();

		//status = freak.PhreeqcRM_ptr->GetSaturation(freak.wsat);
		//freak.temp = freak.PhreeqcRM_ptr->GetTemperature();
		//const std::vector<double> &  tempc = freak.PhreeqcRM_ptr->GetTemperature();
		//if (freak.EK) {status = freak.PhreeqcRM_ptr->GetSpeciesConcentrations(freak.c);}
		//else {		
		//std::vector<double>c(nstot * compn.size());
		//status = freak.PhreeqcRM_ptr->GetConcentrations(c);
		//freak.setC(c);

		std::cout << "c_0 ";for (int i = 0; i < ph_ncomp;i++) { std::cout << freak.c[i * nxyz + 0] << " "; } std::cout << "\n";
		if (gcompn.size() > 0) {
			status = freak.PhreeqcRM_ptr->GetGasCompPressures(freak.g);
			status = freak.PhreeqcRM_ptr->GetGasCompMoles(freak.gm);
			status = freak.PhreeqcRM_ptr->GetGasPhaseVolume(freak.gvol);
			//freak.p = freak.PhreeqcRM_ptr->GetPressure();
			std::cout<<"nb gcomp after "<<gcompn.size()<<" gvol 0 "<<freak.gvol[0]<<" sat "<<freak.wsat[0]<<"\n";
			std::cout << "g_0  ";for (int i = 0;i < ph_gcomp;i++) { std::cout << freak.g[i * nxyz + 0] << " "; } std::cout << "\n";
			std::cout << "gm_0 ";for (int i = 0;i < ph_gcomp;i++) { std::cout << freak.gm[i * nxyz + 0] << " "; } std::cout << "\n";
			std::cout << "p    " << freak.p[0] << " " << "\n";
		}
		status = freak.PhreeqcRM_ptr->GetConcentrations(freak.c);
		std::cout << "c_0 ";for (int i = 0; i < ph_ncomp;i++) { std::cout << freak.c[i * nxyz + 0] << " "; } std::cout << "\n";

		std::cout << "cells run done " << std::endl;
	}
	catch (PhreeqcRMStop)
	{
		std::cerr << "failed with an error in PhreeqcRM." << std::endl;
	#ifdef USE_MPI
		MPI_Abort(MPI_COMM_WORLD, 1);
	#endif
		return IRM_FAIL;
	}
	catch (...)
	{
		std::cerr << "failed with an unhandled exception." << std::endl;
	#ifdef USE_MPI
		MPI_Abort(MPI_COMM_WORLD, 1);
	#endif
		return IRM_FAIL;
	}
	return EXIT_SUCCESS;

};