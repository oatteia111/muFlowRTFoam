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
	PhreeqcRM *PhreeqcRM_ptr;
#ifdef USE_MPI
	MPI_Comm rm_commxx;
#endif
	std::string DB;
	void setDB(std::string db){this->DB = db;}
	std::vector<int> data;
	void setData(std::vector<int> data){this->data = data;}
	std::string ChemFile;
	void setChemFile(std::string chfile){this->ChemFile = chfile;}
	std::vector<double> c;
	void setC(std::vector<double> c){this->c = c;}
	std::vector<double> g;
	void setG(std::vector<double> g){this->g = g;}
	std::vector<double> gm;
	void setGm(std::vector<double> gm){this->gm = gm;}
	std::vector<double> temp;
	void setTemp(std::vector<double> temp){this->temp = temp;}
	std::vector<double> p;
	void setP(std::vector<double> p){this->p = p;}
	std::vector<double> spc;
	void setSpc(std::vector<double> spc){this->spc = spc;}
	//std::vector<double> solu_conc;
	//void setSconc(std::vector<double> solu_conc){this->solu_conc = solu_conc;}
	double tstep=1;
	void setTstep(double tstep){this->tstep = tstep;}
	std::vector<double> poro;
	void setPoro(std::vector<double> poro){this->poro = poro;}
	std::vector<double> wsat;
	void setWsat(std::vector<double> wsat){this->wsat = wsat;}
	std::vector<double> gvol;
	void setGvol(std::vector<double> gvol){this->gvol = gvol;}
	int iGwater=-1;
	void setNwat(int){this->iGwater = iGwater;}
	bool EK;
	int nspc,nselect;
	std::vector<double> diff25;
	std::vector<double> z;
	std::vector<std::string> spcNames;
};
	
int phqInit(my_phq phq_data)
{
	int nxyz=3;int nthreads=2;
	int i,j,ph_ncomp,ph_gcomp,ph_nsolu;
	std::vector<int> ph_init = phq_data.data;
	nxyz=ph_init[0];ph_ncomp=ph_init[1];ph_gcomp=ph_init[2];ph_nsolu=ph_init[3];
	std::cout<<"start of phqInit \n";
	try
	{
	#ifdef USE_MPI
		int argc; char** argv;
		int mpi_myself,mpi_tasks;
		// MPI
		MP_TYPE comm = MPI_COMM_WORLD;
		if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {return EXIT_FAILURE;}
		PhreeqcRM phreeqc_rm(nxyz, MPI_COMM_WORLD);
		phq_data.PhreeqcRM_ptr = &phreeqc_rm;
		if (MPI_Comm_size(MPI_COMM_WORLD, &mpi_tasks) != MPI_SUCCESS)
		{ return EXIT_FAILURE;}

		if (MPI_Comm_rank(MPI_COMM_WORLD, &mpi_myself) != MPI_SUCCESS)
		{
			exit(4);
		}
		if (mpi_myself > 0)
		{
			//phreeqc_rm.SetMpiWorkerCallbackC(worker_tasks_cc);
			//phreeqc_rm.SetMpiWorkerCallbackCookie(&some_data);
			phreeqc_rm.MpiWorker();
			return EXIT_SUCCESS;
		}
	#else
		// OpenMP
		int nthreads = 3;
		PhreeqcRM phreeqc_rm(nxyz, nthreads);
		phq_data.PhreeqcRM_ptr = &phreeqc_rm;
	#endif
		IRM_RESULT status;
		std::cout<<"status reached \n";
		// Set properties
		status = phreeqc_rm.SetErrorHandlerMode(1);
		status = phreeqc_rm.SetComponentH2O(true);
		status = phreeqc_rm.SetRebalanceFraction(0.0);
		status = phreeqc_rm.SetRebalanceByCell(true);
		phreeqc_rm.UseSolutionDensityVolume(false);
		phreeqc_rm.SetPartitionUZSolids(false);
		// Open files
		status = phreeqc_rm.SetFilePrefix("phq1_cpp");
		phreeqc_rm.OpenFiles();
		// Set concentration units
		status = phreeqc_rm.SetUnitsSolution(2);           // 1, mg/L; 2, mol/L; 3, kg/kgs
		status = phreeqc_rm.SetUnitsPPassemblage(1);       // 0, mol/L cell; 1, mol/L water; 2 mol/L rock
		status = phreeqc_rm.SetUnitsExchange(1);           // 0, mol/L cell; 1, mol/L water; 2 mol/L rock
		status = phreeqc_rm.SetUnitsSurface(1);            // 0, mol/L cell; 1, mol/L water; 2 mol/L rock
		status = phreeqc_rm.SetUnitsGasPhase(1);           // 0, mol/L cell; 1, mol/L water; 2 mol/L rock
		status = phreeqc_rm.SetUnitsSSassemblage(1);       // 0, mol/L cell; 1, mol/L water; 2 mol/L rock
		status = phreeqc_rm.SetUnitsKinetics(1);           // 0, mol/L cell; 1, mol/L water; 2 mol/L rock

		std::vector<double> rv(nxyz,1.);
		status = phreeqc_rm.SetRepresentativeVolume(rv);   // RV is one dm3 of medium
		std::cout<<"in phq poro size "<<phq_data.poro.size()<<"\n";
		status = phreeqc_rm.SetPorosity(phq_data.poro);
		std::cout<<"end status \n";
		// Open for writing file
		//std::quoted prf = "phq_out";
		status = phreeqc_rm.SetFilePrefix("phq1");
		phreeqc_rm.OpenFiles();
		//std::cout << "Outfile opened " << "\n";
		//set mask
		std::vector<int> print_chemistry_mask;
		print_chemistry_mask.resize(nxyz, 1);
		//for (int i=0;i<nxyz;i++) {print_chemistry_mask[i]=1;}
		status = phreeqc_rm.SetPrintChemistryMask(print_chemistry_mask);

		//load database and chemistry files and make the first run of the chem file
		status = phreeqc_rm.SetPrintChemistryOn(true, true, true); // workers, initial_phreeqc, utility
		status = phreeqc_rm.LoadDatabase(phq_data.DB);
		std::cout << "dbase opened " << "\n";
		bool workers = true;             // Worker instances do the reaction calculations for transport
		bool initial_phreeqc = true;     // InitialPhreeqc instance accumulates initial and boundary conditions
		bool utility = true;             // Utility instance is available for processing
		status = phreeqc_rm.RunFile(workers, initial_phreeqc,utility, phq_data.ChemFile);

		// Clear contents of workers and utility : seems to create problems with MPI
		/*
		initial_phreeqc = false;
		std::string input = "DELETE; -all";
		status = phreeqc_rm.RunString(workers, initial_phreeqc, utility, input.c_str());
		*/
		int ph_ncomp = phreeqc_rm.FindComponents(); // NECESSARY don't remove
		const std::vector<std::string> &gcompn = phreeqc_rm.GetGasComponents(); //std::cout<<"phq gcomp0 "<<gcomp[0]<<"\n";
		const std::vector<std::string> &compn = phreeqc_rm.GetComponents(); 
		std::cout<<"phq comp "; for (i=0;i<ph_ncomp;i++) {std::cout<<compn[i]<<" ";} std::cout <<"\n";

		const std::vector<int> &f_map = phreeqc_rm.GetForwardMapping();

		// Set array of initial conditions by reading the solution... numbers in data
		std::vector<int> ic1, ic2;
		ic1.resize(nxyz*7, -1);
		//std::vector<double> f1;
		int icnt = 4; // to consider the first 5 numbers : nxyz,ncomp,gcomp,nsolu,units (++ in the first line of the for makes 5 numbers)
		for (i = 0; i < 7; i++) // Solution,Equilibrium phases,Exchange,Surface,Gas phase,Solid solutions,Kinetics
		{
			icnt++;
			int test = phq_data.data[icnt];//std::cout<<"index "<<i<<" " <<test<<"\n";
			if (test == 0)  // if 0 just one value to read
				{icnt++; 
				for (int j=0;j<nxyz;j++) {ic1[i*nxyz+j] = phq_data.data[icnt];
				if (j==1) {std::cout << i << " cell "<< j << " " << phq_data.data[icnt]<< "\n";} }
				}
			else 		
				{for (int j=0;j<nxyz;j++) {icnt++;ic1[i*nxyz+j] = phq_data.data[icnt];
				if (j==1) {std::cout << i << " cell "<< j << " " << phq_data.data[icnt]<< "\n";} }
			}
		}
		//const std::vector<IPhreeqcPhast *> w = phreeqc_rm.GetWorkers();
		//w[0]->AccumulateLine("Delete; -all");
		//int w0 = w[0]->RunAccumulated();

		status = phreeqc_rm.InitialPhreeqc2Module(ic1);

		int ncell = phreeqc_rm.GetChemistryCellCount();
		std::cout << "ncell "<< ncell << " ph_ncomp "<< ph_ncomp<< " starting init phreeqc \n";
		//phreeqc_rm.SetPrintChemistryOn(true,false,false);
		status = phreeqc_rm.SetTemperature(phq_data.temp);
		std::cout << "temp set \n";
		status = phreeqc_rm.RunCells();
		
		std::cout << " cells run \n";
		std::vector<double> temp(nxyz);
		temp = phreeqc_rm.GetTemperature();
		std::cout<<"got temp \n";
		std::vector<double> wsat(nxyz); 
		status = phreeqc_rm.GetSaturation(wsat);
		std::cout<<"got wsat \n";
		std::vector<double> c(nxyz*compn.size());
		status = phreeqc_rm.GetConcentrations(c);
		std::cout<<"got conc  \n";
		
		if (ph_gcomp>0) {
			status = phreeqc_rm.GetSaturation(phq_data.wsat);
			status = phreeqc_rm.GetGasCompPressures(phq_data.g);
			status = phreeqc_rm.GetGasCompMoles(phq_data.gm);
			status = phreeqc_rm.GetGasPhaseVolume(phq_data.gvol);
			std::cout << "in phq poro "<<phq_data.poro[1]<<"\n";
			//std::cout << "c 1-4 ";for (int i=0; i<4;i++){std::cout<< phq_data.c[i*nxyz+1]<<" ";} std::cout<<"\n";
			//std::cout << "c 5-9 ";for (int i=4; i<8;i++){std::cout<< phq_data.c[i*nxyz+1]<<" ";} std::cout<<"\n";
			//std::cout<<"gvol 0 "<<phq_data.gvol[0]<<"gvol 1 "<<phq_data.gvol[1]<<"\n";
			//std::cout<<"g 0 ";for (int i=0;i<ph_gcomp;i++) {std::cout <<phq_data.g[i*nxyz+0]<<" ";} std::cout<<"\n";
			//std::cout<<"g 1 ";for (int i=0;i<ph_gcomp;i++) {std::cout <<phq_data.g[i*nxyz+1]<<" ";} std::cout<<"\n";
			for (size_t i=0;i<gcompn.size();i++) {
				if (gcompn[i]=="H2O(g)") {phq_data.iGwater = i;}
			}
		}
		
	//std::cout << "c ";for (int i=0; i<ph_ncomp;i++){std::cout<< phq_data.c[i*nxyz+0]<<" ";} std::cout<<"\n";
	std::cout <<"end phq init "<<"\n";

		// Get pointer to worker
		int iphreeqc_result;
		const std::vector<IPhreeqcPhast *> w = phreeqc_rm.GetWorkers();
		w[0]->AccumulateLine("Delete; -all");
		iphreeqc_result = w[0]->RunAccumulated();
		// Clean up
		status = phreeqc_rm.CloseFiles();
		status = phreeqc_rm.MpiWorkerBreak();	
	#if defined(USE_MPI)
		MPI_Finalize();
	#endif	
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
	
int phqRun(my_phq &phq_data)
	{
	try
	{
		IRM_RESULT status;
		PhreeqcRM* phreeqc_rm_ptr = phq_data.PhreeqcRM_ptr;
		int nxyz = phq_data.data[0];int ph_gcomp = phq_data.data[2];int ph_nsolu = phq_data.data[3];
		//int nthreads = 3;std::cout<<"threads "<<nthreads<<"\n";
		
		std::vector<int> gcompn;
	status = phreeqc_rm.SetErrorHandlerMode(1);
		int ph_ncomp = phq_data.data[1];//phreeqc_rm.FindComponents(); 
		const std::vector<std::string> &components = phreeqc_rm.GetComponents();
		if (ph_gcomp>0) {const std::vector<std::string> &gcompn = phreeqc_rm.GetGasComponents(); }
		std::cout<<"phqrun compo done";//std::cout<<"phq ncomp sz "<<ph_ncomp<<" gcomp "<<gcompn.size()<<" ws "<<freak.wsat[0]<<"\n";
		std::cout << "in phq nxyz "<<nxyz;//<<" size c "<<freak.c.size()<<" poro "<<freak.poro[1]<<"\n";
		status = phreeqc_rm.SetTemperature(phq_data.temp);std::cout<<"temp done ";
		std::cout << "c_0 ";for (int i=0; i<ph_ncomp;i++){std::cout<< phq_data.c[i*nxyz+2]<<" ";} std::cout<<"\n";
		if (gcompn.size()>0) {
			std::cout<<"nb gcomp before "<<gcompn.size();//<<" gvol 3 "<<freak.gvol[0]<<" sat "<<freak.wsat[0]<<"\n";
			//std::cout<<"g 0 ";for (int i=0;i<ph_gcomp;i++) {std::cout <<freak.g[i*nxyz+0]<<" ";} std::cout<<"\n";
			//std::cout<<"gm 0 ";for (int i=0;i<ph_gcomp;i++) {std::cout <<freak.gm[i*nxyz+0]<<" ";} std::cout<<"\n";
			//status = phreeqc_rm_ptr->SetSaturation(freak.wsat); // keep this order to set values
			status = phreeqc_rm.SetGasPhaseVolume(phq_data.gvol);
			status = phreeqc_rm.SetGasCompMoles(phq_data.gm);
			//status = phreeqc_rm_ptr->SetPressure(freak.p); // does not seem to be usefull 
			}
		else 
			{std::cout<<"phq wsat "<<phq_data.wsat[0];
			//status = phreeqc_rm_ptr->SetSaturation(freak.wsat);
			}
		std::cout<<"wsat done ";
		if (phq_data.EK) {status = phreeqc_rm.SpeciesConcentrations2Module(phq_data.spc);}
		else {status = phreeqc_rm.SetConcentrations(phq_data.c);}
		std::cout<<"c done ";
		status = phreeqc_rm.SetTimeStep(phq_data.tstep);std::cout<<"tstep "<<phq_data.tstep<<"\n";
		//----- run here ------
		status = phreeqc_rm.RunCells();
		
		//status = phreeqc_rm_ptr->GetSaturation(freak.wsat);
		//freak.temp = phreeqc_rm_ptr->GetTemperature();
		if (phq_data.EK) {status = phreeqc_rm.GetSpeciesConcentrations(phq_data.c);}
		else {status = phreeqc_rm.GetConcentrations(phq_data.c);}
	std::cout << "c 0 ";for (int i=0; i<ph_ncomp;i++){std::cout<< phq_data.c[i*nxyz]<<" ";} std::cout<<"\n";
		if (gcompn.size()>0) {
			status = phreeqc_rm.GetGasCompPressures(phq_data.g);
			status = phreeqc_rm.GetGasCompMoles(phq_data.gm);
			status = phreeqc_rm.GetGasPhaseVolume(phq_data.gvol);
			phq_data.p = phreeqc_rm.GetPressure();
			}
		//writes to the dump file
		//bool dump_on = true;
		//bool append = false;
		//status = phreeqc_rm_ptr->SetDumpFileName("phreeqc.dmp");
		//status = phreeqc_rm_ptr->DumpModule(dump_on, append);
	//phreeqc_rm_ptr->SetSelectedOutputOn(true);
	std::cout << "c 0 ";for (int i=0; i<ph_ncomp;i++){std::cout<< phq_data.c[i*nxyz+0]<<" ";} std::cout<<"\n";
	if (ph_gcomp>0) {
	//std::cout<<"nb gcomp after "<<gcompn.size()<<" gvol 0 "<<freak.gvol[0]<<" sat "<<freak.wsat[0]<<"\n";
		std::cout<<"g 0 ";for (int i=0;i<ph_gcomp;i++) {std::cout <<phq_data.g[i*nxyz+0]<<" ";} std::cout<<"\n";
		std::cout<<"gm 0 ";for (int i=0;i<ph_gcomp;i++) {std::cout <<phq_data.gm[i*nxyz+0]<<" ";} std::cout<<"\n";
		std::cout<<"p ";for (int i=0;i<5;i++) {std::cout<<phq_data.p[i]<<" ";} std::cout<<"\n";
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

int i,j,nxyz,ph_ncomp,ph_gcomp,ph_nsolu;



int main()
{
	my_phq phq_data;
	phq_data.DB="phreeqc.dat";
	phq_data.ChemFile = "initChem.pqi";
	char sep ='/';std::cout<<"start \n";
	std::ifstream inputInit{cur_dir+sep+"phqinit.txt"};
	std::vector<int> ph_init{std::istream_iterator<int>{inputInit}, {}}; 
	phq_data.setData(ph_init);
	nxyz=phq_data.data[0];std::cout<<"nxyz "<<nxyz<<"\n";
	std::vector<double> rv(nxyz, 1.);
	std::vector<double> poro(nxyz, .25);for (i=0;i<nxyz;i++) {poro[i]=0.25;}
	std::vector<double> temp(nxyz, 25);
	//std::vector<int> data={nxyz,ph_ncomp,0,ph_nsolu,1, 0,-1, 0,-1,0,-1,0,-1,0,-1,0,-1,0,-1};
	phq_data.setTemp(temp);//std::cout("temp size",temp.size());
	phq_data.setPoro(poro);std::cout<<"poro size "<<poro.size()<<"\n";
	int a0=phqInit(phq_data);
	//phq_data.setC(c);
	std::cout<<"init done\n";
	a0 = phqRun(phq_data.PhreeqcRM_ptr);
}

