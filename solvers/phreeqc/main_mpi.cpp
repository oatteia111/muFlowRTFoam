
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

#include <unistd.h>
#define GetCurrentDir getcwd

std::string get_current_dir() {
   char buff[FILENAME_MAX]; //create string buffer to hold path
   GetCurrentDir( buff, FILENAME_MAX );
   std::string current_working_dir(buff);
   return current_working_dir;
}
std::string cur_dir = get_current_dir();

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
	int iGwater = -1;
	void setNwat(int) { this->iGwater = iGwater; }
	bool EK;
	int nspc, nselect;
	std::vector<double> diff25;
	std::vector<double> z;
	std::vector<std::string> spcNames;
};

int phqInit(my_phq & phq_data)
{
	int nxyz = 3;int nthreads = 2;
	int i, j, ph_ncomp, ph_gcomp, ph_nsolu;
	std::vector<int> ph_init = phq_data.data;
	nxyz = ph_init[0];ph_ncomp = ph_init[1];ph_gcomp = ph_init[2];ph_nsolu = ph_init[3];
	PhreeqcRM* phreeqc_rm_ptr = phq_data.PhreeqcRM_ptr;
	std::cout << "start of phqInit " << nxyz << std::endl;
	try
	{
		IRM_RESULT status;
		// Set properties
		status = phreeqc_rm_ptr->SetErrorHandlerMode(1);
		status = phreeqc_rm_ptr->SetComponentH2O(true);
		status = phreeqc_rm_ptr->SetRebalanceFraction(0.0);
		status = phreeqc_rm_ptr->SetRebalanceByCell(true);
		phreeqc_rm_ptr->UseSolutionDensityVolume(false);
		phreeqc_rm_ptr->SetPartitionUZSolids(false);
		std::cerr << "before  SetFilePrefix " << std::endl;
		// Open files
		status = phreeqc_rm_ptr->SetFilePrefix("phq1_cpp");
		phreeqc_rm_ptr->OpenFiles();
		// Set concentration units
		status = phreeqc_rm_ptr->SetUnitsSolution(2);           // 1, mg/L; 2, mol/L; 3, kg/kgs
		status = phreeqc_rm_ptr->SetUnitsPPassemblage(1);       // 0, mol/L cell; 1, mol/L water; 2 mol/L rock
		status = phreeqc_rm_ptr->SetUnitsExchange(1);           // 0, mol/L cell; 1, mol/L water; 2 mol/L rock
		status = phreeqc_rm_ptr->SetUnitsSurface(1);            // 0, mol/L cell; 1, mol/L water; 2 mol/L rock
		status = phreeqc_rm_ptr->SetUnitsGasPhase(1);           // 0, mol/L cell; 1, mol/L water; 2 mol/L rock
		status = phreeqc_rm_ptr->SetUnitsSSassemblage(1);       // 0, mol/L cell; 1, mol/L water; 2 mol/L rock
		status = phreeqc_rm_ptr->SetUnitsKinetics(1);           // 0, mol/L cell; 1, mol/L water; 2 mol/L rock

		status = phreeqc_rm_ptr->SetRepresentativeVolume(phq_data.rv);   // RV is one dm3 of medium
		std::cerr << "before  SetPorosity " << phq_data.poro.size() << std::endl;
		status = phreeqc_rm_ptr->SetPorosity(phq_data.poro);
		std::cout << "end status \n";
		// Open for writing file
		status = phreeqc_rm_ptr->SetFilePrefix("phq1");
		phreeqc_rm_ptr->OpenFiles();
		//set mask
		std::vector<int> print_chemistry_mask;
		print_chemistry_mask.resize(nxyz, 1);
		status = phreeqc_rm_ptr->SetPrintChemistryMask(print_chemistry_mask);

		//load database and chemistry files and make the first run of the chem file
		status = phreeqc_rm_ptr->SetPrintChemistryOn(true, true, true); // workers, initial_phreeqc, utility
		status = phreeqc_rm_ptr->LoadDatabase(phq_data.DB);
		std::cout << "dbase opened " << "\n";
		bool workers = true;             // Worker instances do the reaction calculations for transport
		bool initial_phreeqc = true;     // InitialPhreeqc instance accumulates initial and boundary conditions
		bool utility = true;             // Utility instance is available for processing
		status = phreeqc_rm_ptr->RunFile(workers, initial_phreeqc, utility, phq_data.ChemFile);
		std::cerr << "After RunFile." << std::endl;

		int ph_ncomp = phreeqc_rm_ptr->FindComponents(); // NECESSARY don't remove
		const std::vector<std::string>& gcompn = phreeqc_rm_ptr->GetGasComponents(); //std::cout<<"phq gcomp0 "<<gcomp[0]<<"\n";
		const std::vector<std::string>& compn = phreeqc_rm_ptr->GetComponents();
		std::cout << "phq comp "; for (i = 0;i < ph_ncomp;i++) { std::cout << compn[i] << " "; } std::cout << "\n";

		const std::vector<int>& f_map = phreeqc_rm_ptr->GetForwardMapping();

		// Set array of initial conditions by reading the solution... numbers in data
		std::vector<int> ic1, ic2;
		ic1.resize(nxyz * 7, -1);
		int icnt = 4; // to consider the first 5 numbers : nxyz,ncomp,gcomp,nsolu,units (++ in the first line of the for makes 5 numbers)
		for (i = 0; i < 7; i++)  //Solution,Equilibrium phases,Exchange,Surface,Gas phase,Solid solutions,Kinetics
		{
			icnt++;
			int test = phq_data.data[icnt];//std::cout<<"index "<<i<<" " <<test<<"\n";
			if (test == 0)  // if 0 just one value to read
			{
				icnt++;
				for (int j = 0;j < nxyz;j++) {
					ic1[i * nxyz + j] = phq_data.data[icnt];
					//if (j == 1) { std::cout << i << " cell " << j << " " << phq_data.data[icnt] << "\n"; }
				}
			}
			else
			{
				for (int j = 0;j < nxyz;j++) {
					icnt++;ic1[i * nxyz + j] = phq_data.data[icnt];
					//if (j == 1) { std::cout << i << " cell " << j << " " << phq_data.data[icnt] << "\n"; }
				}
			}
		}
		status = phreeqc_rm_ptr->InitialPhreeqc2Module(ic1);

		int ncell = phreeqc_rm_ptr->GetChemistryCellCount();
		std::cout << "ncell " << ncell << " ph_ncomp " << ph_ncomp << " starting init phreeqc \n";
		//phreeqc_rm_ptr->SetPrintChemistryOn(true,false,false);
		status = phreeqc_rm_ptr->SetTemperature(phq_data.temp);
		status = phreeqc_rm_ptr->RunCells();
		std::cout << " cells run \n";
		std::vector<double> temp(nxyz);
		temp = phreeqc_rm_ptr->GetTemperature();
		std::vector<double> wsat(nxyz);
		status = phreeqc_rm_ptr->GetSaturation(wsat);
		std::vector<double> c(nxyz * compn.size());
		status = phreeqc_rm_ptr->GetConcentrations(c);
		phq_data.setC(c);
		std::cerr << "After GetConcentrations. " << std::endl;

		if (ph_gcomp > 0) {
			status = phreeqc_rm_ptr->GetSaturation(phq_data.wsat);
			status = phreeqc_rm_ptr->GetGasCompPressures(phq_data.g);
			status = phreeqc_rm_ptr->GetGasCompMoles(phq_data.gm);
			status = phreeqc_rm_ptr->GetGasPhaseVolume(phq_data.gvol);
			std::cout << "in phq poro " << phq_data.poro[1] << "\n";
			assert(compn.size() == (size_t)phreeqc_rm_ptr->GetComponentCount());
			for (size_t i = 0; i < compn.size(); i++)
			{
				std::cerr << compn[i] << " " << phq_data.c[i * nxyz + 1] << std::endl;
			}
			//std::cout<<"gvol 0 "<<phq_data.gvol[0]<<"gvol 1 "<<phq_data.gvol[1]<<"\n";
			//std::cout<<"g 0 ";for (int i=0;i<ph_gcomp;i++) {std::cout <<phq_data.g[i*nxyz+0]<<" ";} std::cout<<"\n";
			//std::cout<<"g 1 ";for (int i=0;i<ph_gcomp;i++) {std::cout <<phq_data.g[i*nxyz+1]<<" ";} std::cout<<"\n";
			for (size_t i = 0;i < gcompn.size();i++) {
				if (gcompn[i] == "H2O(g)") { phq_data.iGwater = i; }
			}
		}
		//std::cout << "c ";for (int i=0; i<ph_ncomp;i++){std::cout<< phq_data.c[i*nxyz+0]<<" ";} std::cout<<"\n";
		std::cout << "end phq init " << "\n";

		// Get pointer to worker
		int iphreeqc_result;
		const std::vector<IPhreeqcPhast*> w = phreeqc_rm_ptr->GetWorkers();
		w[0]->AccumulateLine("Delete; -all");
		iphreeqc_result = w[0]->RunAccumulated();
		status = phreeqc_rm_ptr->CloseFiles();

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
	std::cerr << "End phqInit" << std::endl;
	return EXIT_SUCCESS;
};

	//--------------------------- RUN ------------------------------
	
int phqRun(my_phq & freak)
	{
	try
	{
		IRM_RESULT status;
		PhreeqcRM* phreeqc_rm_ptr = freak.PhreeqcRM_ptr;
		
		int ph_ncomp = freak.PhreeqcRM_ptr->FindComponents(); 
		const std::vector<std::string> &components = phreeqc_rm_ptr->GetComponents();
		const std::vector<std::string> &gcompn = phreeqc_rm_ptr->GetGasComponents(); std::cout<<"phq ncomp sz "<<ph_ncomp<<" gcomp "<<gcompn.size()<<" ws "<<freak.wsat[0]<<"\n";
		int nxyz = freak.data[0];int ph_gcomp = freak.data[2];int ph_nsolu = freak.data[3];
		std::cout << "in phq nxyz "<<nxyz<<" size c "<<freak.c.size()<<" poro "<<freak.poro[2]<<"\n";
		status = phreeqc_rm_ptr->SetTemperature(freak.temp);std::cout<<"temp done ";
		std::cout << "c_0 ";for (int i=0; i<ph_ncomp;i++){std::cout<< freak.c[i*nxyz+2]<<" ";} std::cout<<"\n";
		if (gcompn.size()>0) {
			std::cout<<"nb gcomp before "<<gcompn.size();//<<" gvol 3 "<<freak.gvol[0]<<" sat "<<freak.wsat[0]<<"\n";
			//std::cout<<"g 0 ";for (int i=0;i<ph_gcomp;i++) {std::cout <<freak.g[i*nxyz+0]<<" ";} std::cout<<"\n";
			//std::cout<<"gm 0 ";for (int i=0;i<ph_gcomp;i++) {std::cout <<freak.gm[i*nxyz+0]<<" ";} std::cout<<"\n";
			//status = phreeqc_rm_ptr->SetSaturation(freak.wsat); // keep this order to set values
			status = phreeqc_rm_ptr->SetGasPhaseVolume(freak.gvol);
			status = phreeqc_rm_ptr->SetGasCompMoles(freak.gm);
			//status = phreeqc_rm_ptr->SetPressure(freak.p); // does not seem to be usefull 
			}
		else 
			{std::cout<<"phq wsat "<<freak.wsat[0];
			status = phreeqc_rm_ptr->SetSaturation(freak.wsat);
			}
		std::cout<<"wsat done ";
		if (freak.EK) {status = phreeqc_rm_ptr->SpeciesConcentrations2Module(freak.spc);}
		else {		status = phreeqc_rm_ptr->SetConcentrations(freak.c);}
		std::cout<<"c done ";
		status = phreeqc_rm_ptr->SetTimeStep(freak.tstep);std::cout<<"tstep "<<freak.tstep<<"\n";
		//----- run here ------
		status = phreeqc_rm_ptr->RunCells();
		
		status = phreeqc_rm_ptr->GetSaturation(freak.wsat);
		//freak.temp = phreeqc_rm_ptr->GetTemperature();
		//const std::vector<double> &  tempc = phreeqc_rm_ptr->GetTemperature();

		if (freak.EK) {status = phreeqc_rm_ptr->GetSpeciesConcentrations(freak.c);}
		else {		status = phreeqc_rm_ptr->GetConcentrations(freak.c);}
	std::cout << "c 0 ";for (int i=0; i<ph_ncomp;i++){std::cout<< freak.c[i*nxyz]<<" ";} std::cout<<"\n";
		if (gcompn.size()>0) {
			status = phreeqc_rm_ptr->GetGasCompPressures(freak.g);
			status = phreeqc_rm_ptr->GetGasCompMoles(freak.gm);
			status = phreeqc_rm_ptr->GetGasPhaseVolume(freak.gvol);
			//freak.p = phreeqc_rm_ptr->GetPressure();
			}
		//writes to the dump file
		//bool dump_on = true;
		//bool append = false;
		//status = phreeqc_rm_ptr->SetDumpFileName("phreeqc.dmp");
		//status = phreeqc_rm_ptr->DumpModule(dump_on, append);
	//phreeqc_rm_ptr->SetSelectedOutputOn(true);
	std::cout << "c 0 ";for (int i=0; i<ph_ncomp;i++){std::cout<< freak.c[i*nxyz+0]<<" ";} std::cout<<"\n";
	if (ph_gcomp>0) {
	//std::cout<<"nb gcomp after "<<gcompn.size()<<" gvol 0 "<<freak.gvol[0]<<" sat "<<freak.wsat[0]<<"\n";
		std::cout<<"g 0 ";for (int i=0;i<ph_gcomp;i++) {std::cout <<freak.g[i*nxyz+0]<<" ";} std::cout<<"\n";
		std::cout<<"gm 0 ";for (int i=0;i<ph_gcomp;i++) {std::cout <<freak.gm[i*nxyz+0]<<" ";} std::cout<<"\n";
		std::cout<<"p ";for (int i=0;i<5;i++) {std::cout<<freak.p[i]<<" ";} std::cout<<"\n";
		}
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

int i, j, nxyz, ph_ncomp, ph_gcomp, ph_nsolu;

int main(int argc, char* argv[])
{
	my_phq freak;
	cur_dir = get_current_dir();
	//freak.EK=false;
	char sep = '/';std::cout << "start \n";
	std::ifstream inputInit{cur_dir + sep + "phqinit.txt"};
	std::vector<int> ph_init{std::istream_iterator<int>{inputInit}, {}};
	nxyz = ph_init[0];ph_ncomp=ph_init[1];ph_gcomp=ph_init[2];ph_nsolu=ph_init[3];
	int mpi_myself, mpi_tasks;
#ifdef USE_MPI
	// MPI
	MP_TYPE comm = MPI_COMM_WORLD;
	if (MPI_Init(&argc, &argv) != MPI_SUCCESS) { return EXIT_FAILURE; }
	if (MPI_Comm_size(MPI_COMM_WORLD, &mpi_tasks) != MPI_SUCCESS)
	{
		return EXIT_FAILURE;
	}

	if (MPI_Comm_rank(MPI_COMM_WORLD, &mpi_myself) != MPI_SUCCESS)
	{
		exit(4);
	}
	std::cerr << "MPI started. " << mpi_myself << std::endl;
	if (mpi_myself == 0)
	{
		freak.setDB ("phreeqc.dat");
		freak.setChemFile("initChem.pqi");
		freak.setData(ph_init);
	}
	MPI_Bcast(&nxyz, 1, MPI_INT, 0, MPI_COMM_WORLD);
	PhreeqcRM phreeqc_rm(nxyz, MPI_COMM_WORLD);
	freak.PhreeqcRM_ptr = &phreeqc_rm;
	if (mpi_myself > 0)
	{
		std::cerr << "Started MPI worker " << mpi_myself << std::endl;
		phreeqc_rm.MpiWorker();
		MPI_Finalize();
		return EXIT_SUCCESS;
	}
	std::cerr << "Started MPI root " << mpi_myself << std::endl;
#else
	// OpenMP
	int nthreads = 3;
	PhreeqcRM phreeqc_rm(nxyz, nthreads);
	freak.PhreeqcRM_ptr = &phreeqc_rm;
#endif
	// All processing by root from here on
	std::vector<double> rv(nxyz, 1.);freak.setRV(rv);
	std::vector<double> poro(nxyz, .25);freak.setPoro(poro);
	std::vector<double> wsat(nxyz, 0.9);freak.setWsat(wsat);
	std::vector<double> temp(nxyz, 25.);freak.setTemp(temp);
	std::vector<double> p(nxyz, 1.);freak.setP(p);
	std::vector<double> g(nxyz*2, 1.);freak.setP(g);
	std::vector<double> gm(nxyz*2, 1.);freak.setP(gm);
	std::vector<double> gvol(nxyz, 0.1);freak.setGvol(gvol);

	int a0 = phqInit(freak);
	//std::vector<double> c(nxyz*ph_nsolu);
	//c=freak.c;
	std::cout << "init done c0 "<<freak.c[0]<<"\n";
	a0 = phqRun(freak);

#ifdef USE_MPI
	phreeqc_rm.MpiWorkerBreak();
	std::cerr << "After MpiWorkerBreak." << std::endl;
	int status = MPI_Finalize();
	std::cerr << "After MPI_Finalize root. " << status << std::endl;
#endif

	exit(EXIT_SUCCESS);
}

