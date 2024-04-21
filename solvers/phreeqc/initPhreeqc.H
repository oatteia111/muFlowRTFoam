#if defined(USE_MPI)
#include <mpi.h>
#endif
#include "PhreeqcRM.h"
#include "IPhreeqc.hpp"
#include "IPhreeqcPhast.h"
//int nthreads = 3;

int nxyz,ph_ncomp,ph_gcomp,ph_nsolu;
double my_basic_callback(double x1, double x2, const char *str, void *cookie);
void register_basic_callback(void *cookie);

////////////////////// start phreeqc
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
	std::vector<double> bc_conc;
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
	std::vector<double> foc;
	void setFoc(std::vector<double> foc){this->foc = foc;}
	std::vector<std::string> comp;
	std::vector<std::string> gcomp;
	int iGwater=-1;
	void setNwat(int){this->iGwater = iGwater;}
	bool EK;
	int nspc,nselect;
	std::vector<double> diff25;
	std::vector<double> z;
	std::vector<std::string> spcNames;

	//std::vector<std::string> complist;
    //void setCompList(std::vector<std::string>) {this->complist = complist;}
	//this->complist = complist;
};
int phqInit(my_phq & freak)
{
	nxyz = freak.data[0];ph_ncomp = freak.data[1];ph_gcomp = freak.data[2];
	ph_nsolu=freak.data[3];std::cerr << " ph nxyz "<< nxyz <<"\n";
	int i,j;
	if (ph_ncomp<0) {bool EK = true;freak.EK = EK;}
	else {std::vector<double> c(nxyz*ph_ncomp,0.);freak.c = c;}
	
	#ifdef USE_MPI
		PhreeqcRM* phreeqc_rm = freak.PhreeqcRM_ptr;
	#else
		
		int nthreads;
		PhreeqcRM *phreeqc_rm = new PhreeqcRM (nxyz, nthreads);
		freak.PhreeqcRM_ptr = phreeqc_rm;		
	#endif
	
	try
	{	
		IRM_RESULT status;
		status = phreeqc_rm->SetErrorHandlerMode(1);
		status = phreeqc_rm->SetScreenOn(false);
		status = phreeqc_rm->SetComponentH2O(true); // to get water,excess H an dexcess O transported
		status = phreeqc_rm->SetRebalanceFraction(0.0);
		status = phreeqc_rm->SetRebalanceByCell(true);
		phreeqc_rm->UseSolutionDensityVolume(false);
		phreeqc_rm->SetPartitionUZSolids(false);
		status = phreeqc_rm->SetUnitsSolution(2);           // 1, mg/L; 2, mol/L; 3, kg/kgs
		status = phreeqc_rm->SetUnitsPPassemblage(1);       // 0, mol/L cell; 1, mol/L water; 2 mol/L rock
		status = phreeqc_rm->SetUnitsExchange(1);           // 0, mol/L cell; 1, mol/L water; 2 mol/L rock
		status = phreeqc_rm->SetUnitsSurface(1);            // 0, mol/L RV; 1, mol/L water; 2 mol/L rock
		status = phreeqc_rm->SetUnitsGasPhase(0);           // 0, mol/L RV; 1, mol/L water; 2 mol/L rock
		status = phreeqc_rm->SetUnitsSSassemblage(1);       // 0, mol/L cell; 1, mol/L water; 2 mol/L rock
		status = phreeqc_rm->SetUnitsKinetics(0);           // 0, mol/L cell; 1, mol/L water; 2 mol/L rock
		if (freak.EK) {status = phreeqc_rm->SetSpeciesSaveOn(true);}
				
		std::vector<double> rv(nxyz, 1.);
		status = phreeqc_rm->SetRepresentativeVolume(rv);   // RV is one dm3 of medium
		std::cout<<"poro size "<<freak.poro.size()<<"\n";
		status = phreeqc_rm->SetPorosity(freak.poro);

		// Open for writing file
		//std::quoted prf = "phq_out";
		status = phreeqc_rm->SetFilePrefix("q");
		phreeqc_rm->OpenFiles();
		//std::cout << "Outfile opened " << "\n";
		//set mask
		std::vector<int> print_chemistry_mask;
		print_chemistry_mask.resize(nxyz, 1);
		//for (int i=0;i<nxyz;i++) {print_chemistry_mask[i]=1;}
		status = phreeqc_rm->SetPrintChemistryMask(print_chemistry_mask);

		//load database and chemistry files and make the first run of the chem file
		status = phreeqc_rm->SetPrintChemistryOn(true, true, true); // workers, initial_phreeqc, utility
		status = phreeqc_rm->LoadDatabase(freak.DB);
		std::cout << "dbase opened " << "\n";
		bool workers = true;             // Worker instances do the reaction calculations for transport
		bool initial_phreeqc = true;     // InitialPhreeqc instance accumulates initial and boundary conditions
		bool utility = true;             // Utility instance is available for processing
		status = phreeqc_rm->RunFile(workers, initial_phreeqc,utility, freak.ChemFile);

		// Clear contents of workers and utility
		initial_phreeqc = false;
		std::string input = "DELETE; -all";
		status = phreeqc_rm->RunString(workers, initial_phreeqc, utility, input.c_str());
		
		int ph_ncomp = phreeqc_rm->FindComponents(); // NECESSARY don't remove
		const std::vector<std::string> &gcomp = phreeqc_rm->GetGasComponents(); //std::cout<<"phq gcomp0 "<<gcomp[0]<<"\n";
		const std::vector<std::string> &comp = phreeqc_rm->GetComponents(); 
		freak.comp = comp;freak.gcomp = gcomp;
		std::cout<<"phq comp "; for (i=0;i<ph_ncomp;i++) {std::cout<<comp[i]<<" ";} std::cout <<"\n";

		// Set array of initial conditions by reading the solution... numbers in data
		std::vector<int> ic1;
		ic1.resize(nxyz*7, -1);std::cout<<"ic1 size "<<ic1.size()<<"\n";
		//std::vector<double> f1;
		int icnt = 4; // to consider the first 5 numbers : nxyz,ncomp,gcomp,nsolu,units (++ in the first line of the for makes 5 numbers)
		for (i = 0; i < 7; i++) // Solution,Equilibrium phases,Exchange,Surface,Gas phase,Solid solutions,Kinetics
		{
			icnt++;
			int test = freak.data[icnt];//std::cout<<"index "<<i<<" " <<test<<"\n";
			if (test == 0)  // if 0 just one value to read
				{icnt++; 
				for (int j=0;j<nxyz;j++) {ic1[i*nxyz+j] = freak.data[icnt];
				if (j==10) {std::cout << i << " cell "<< j << " " << freak.data[icnt]<< "\n";} }
				}
			else 		
				{for (int j=0;j<nxyz;j++) {icnt++;ic1[i*nxyz+j] = freak.data[icnt];
				if (j==10) {std::cout << i << " cell "<< j << " " << freak.data[icnt]<< "\n";} }
			}
		}
		std::cout << "end fill \n";
		/*
		const std::vector<IPhreeqcPhast *> w = phreeqc_rm->GetWorkers();
		w[0]->AccumulateLine("Delete; -all");
		int w0 = w[0]->RunAccumulated();
		*/
		//set initial conditions
		status = phreeqc_rm->InitialPhreeqc2Module(ic1);
		std::cout << "end init Ic \n";
		
		// Demonstrate add to Basic: Set a function for Basic CALLBACK after LoadDatabase
		register_basic_callback(&freak);

		// run initial conditions
		int ncell = phreeqc_rm->GetChemistryCellCount();
		std::cout << "ncell "<< ncell << " ph_ncomp "<< ph_ncomp<< " starting init phreeqc \n";
		status = phreeqc_rm->SetTemperature(freak.temp);std::cout<<"temp "<<freak.temp[0]<<"\n";
		
		status = phreeqc_rm->RunCells();
		
		std::cout << " cells run \n";
		status = phreeqc_rm->GetConcentrations(freak.c);
		
		//set boundary condition
		//std::vector<double> bc_conc;
		std::vector<int> bc1;
		bc1.resize(ph_nsolu,0);
		for (i=0;i<ph_nsolu;i++) {bc1[i]=i;}
		status = phreeqc_rm->InitialPhreeqc2Concentrations(freak.bc_conc, bc1);
	
		if (freak.EK) {
			int nspc = phreeqc_rm->GetSpeciesCount(); //std::cout<<"phq eK nspc "<<EK<<" "<<nspc<<"\n";
			freak.nspc = nspc;
			status = phreeqc_rm->GetSpeciesConcentrations(freak.spc);
			const std::vector <double> & diff25 = phreeqc_rm->GetSpeciesD25();
			freak.diff25 = diff25;
			const std::vector <double> & z = phreeqc_rm->GetSpeciesZ();
			freak.z = z;
			const std::vector<std::string> & spcNames = phreeqc_rm->GetSpeciesNames();
			freak.spcNames = spcNames;
		}
		if (ph_gcomp>0) {
			status = phreeqc_rm->GetSaturation(freak.wsat);
			status = phreeqc_rm->GetGasCompPressures(freak.g);
			status = phreeqc_rm->GetGasCompMoles(freak.gm);
			status = phreeqc_rm->GetGasPhaseVolume(freak.gvol);
			std::cout << "in phq poro "<<freak.poro[1]<<"\n";
			std::cout << "c 1-4 ";for (int i=0; i<4;i++){std::cout<< freak.c[i*nxyz+1]<<" ";} std::cout<<"\n";
			std::cout << "c 5-9 ";for (int i=4; i<8;i++){std::cout<< freak.c[i*nxyz+1]<<" ";} std::cout<<"\n";
			//std::cout<<"Vm "<<freak.poro[0]/freak.gm[0]<<" "<<freak.poro[0]/freak.gm[1]<<"\n";
			//std::cout<<"gvol 0 "<<freak.gvol[0]<<"gvol 1 "<<freak.gvol[1]<<"\n";
			//std::cout<<"g 0 ";for (int i=0;i<ph_gcomp;i++) {std::cout <<freak.g[i*nxyz+0]<<" ";} std::cout<<"\n";
			//std::cout<<"g 1 ";for (int i=0;i<ph_gcomp;i++) {std::cout <<freak.g[i*nxyz+1]<<" ";} std::cout<<"\n";
			for (i=0;i<gcomp.size();i++) {
				if (gcomp[i]=="H2O(g)") {freak.iGwater = i;}
			}
		}
	std::cout << "c ";for (int i=0; i<ph_ncomp;i++){std::cout<< freak.c[i*nxyz+0]<<" ";} std::cout<<"\n";
	std::cout <<"end phq init "<<"\n";

		/* Set boundary condition
		std::vector<double> bc_conc;
		std::vector<double> bc_f1; //bc_conc, 
		std::vector<int> bc1, bc2;
		int nbound = 2;
		bc1.resize(nbound, 0);                      // solution 0 from Initial IPhreeqc instance to boundary condition
		bc2.resize(nbound, 1);                     // no bc2 solution for mixing
		bc_f1.resize(nbound, 1.0);                  // mixing fraction for bc1
		status = phreeqc_rm->InitialPhreeqc2Concentrations(freak.solu_conc, bc1, bc2, bc_f1); */
		
		//printing to file
		status = phreeqc_rm->CloseFiles();
		//status = phreeqc_rm->SetPrintChemistryOn(true, false, true); // workers, initial_phreeqc, utility
	
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
	
	//--------------------------- RUN ------------------------------

int phqRun(my_phq & freak)
{
	try
	{
	IRM_RESULT status;
	//PhreeqcRM* phreeqc_rm = freak.PhreeqcRM_ptr;
	
	status = freak.PhreeqcRM_ptr->SetErrorHandlerMode(1);
	//status = freak.PhreeqcRM_ptr->SetFilePrefix("phq_run");
	//freak.PhreeqcRM_ptr->OpenFiles();
	//std::vector<int> print_chemistry_mask;
	//print_chemistry_mask.resize(nxyz, 1);//print_chemistry_mask[0]=1;print_chemistry_mask[50]=1;
	//status = freak.PhreeqcRM_ptr->SetPrintChemistryMask(print_chemistry_mask); 
	//status = freak.PhreeqcRM_ptr->SetPrintChemistryOn(true, false, false); // workers, initial_phreeqc, utility
	const std::vector<std::string> &components = freak.PhreeqcRM_ptr->GetComponents();
	const std::vector<std::string> &gcomp = freak.PhreeqcRM_ptr->GetGasComponents(); std::cout<<"phq ngcomp sz "<<gcomp.size()<<" phq ws "<<freak.wsat[0]<<"\n";
	nxyz = freak.data[0];
	std::cout << "in phq nxyz "<<nxyz<<" size c "<<freak.c.size()<<" poro "<<freak.poro[2]<<"\n";
	
	status = freak.PhreeqcRM_ptr->SetTemperature(freak.temp);
	std::cout<<"temp 0-5 ";for (int j=0; j<5;j++){std::cout<<freak.temp[j]<<" ";} std::cout<<"\n";
	std::cout<<"press 0-5 ";for (int j=0; j<5;j++){std::cout<<freak.p[j]<<" ";} std::cout<<"\n";
	//for (int j=0;j<nxyz;j++) {float c0=freak.c[6*nxyz+j];if ((c0>4e-4)&&(c0<6e-4)) {std::cout << "cp "<<j<<" ";for (int i=0; i<ph_ncomp;i++){std::cout<< freak.c[i*nxyz+j]<<" ";} std::cout<<"\n";} }
	if (gcomp.size()>0) {
		std::cout<<"nb gcomp before "<<gcomp.size()<<"\n";//<<" gvol 3 "<<freak.gvol[0]<<" sat "<<freak.wsat[0]<<"\n";
		//std::cout<<"g 0 ";for (int i=0;i<ph_gcomp;i++) {std::cout <<freak.g[i*nxyz+0]<<" ";} std::cout<<"\n";
		//std::cout<<"gm 0 ";for (int i=0;i<ph_gcomp;i++) {std::cout <<freak.gm[i*nxyz+0]<<" ";} std::cout<<"\n";
		status = freak.PhreeqcRM_ptr->SetSaturation(freak.wsat); // keep this order to set values
		status = freak.PhreeqcRM_ptr->SetGasPhaseVolume(freak.gvol);
		status = freak.PhreeqcRM_ptr->SetGasCompMoles(freak.gm);
		//status = freak.PhreeqcRM_ptr->SetPressure(freak.p); // does not seem to be usefull but we left it
		}
	else 
		{std::cout<<"phq wsat "<<freak.wsat[0];
		//status = freak.PhreeqcRM_ptr->SetSaturation(freak.wsat);
		}
	std::cout<<"wsat done ";
	if (freak.EK) {status = freak.PhreeqcRM_ptr->SpeciesConcentrations2Module(freak.spc);}
	else {status = freak.PhreeqcRM_ptr->SetConcentrations(freak.c);}
	std::cout<<"c done ";
	status = freak.PhreeqcRM_ptr->SetTimeStep(freak.tstep);std::cout<<"tstep "<<freak.tstep<<"\n";
	//----- run here ------
	status = freak.PhreeqcRM_ptr->RunCells();
	
	//status = freak.PhreeqcRM_ptr->GetSaturation(freak.wsat);
	//status = freak.PhreeqcRM_ptr->GetTemperature(freak.temp);
	if (freak.EK) {status = freak.PhreeqcRM_ptr->GetSpeciesConcentrations(freak.c);}
	else {status = freak.PhreeqcRM_ptr->GetConcentrations(freak.c);}
	std::cout << "c 0 ";for (int i=0; i<ph_ncomp;i++){std::cout<< freak.c[i*nxyz]<<" ";} std::cout<<"\n";
	std::cout << "c end ";for (int i=0; i<ph_ncomp;i++){std::cout<< freak.c[i*nxyz+nxyz-1]<<" ";} std::cout<<"\n";
	if (gcomp.size()>0) {
		status = freak.PhreeqcRM_ptr->GetGasCompPressures(freak.g);
		status = freak.PhreeqcRM_ptr->GetGasCompMoles(freak.gm);
		status = freak.PhreeqcRM_ptr->GetGasPhaseVolume(freak.gvol);
		freak.p = freak.PhreeqcRM_ptr->GetPressure();
		}
	//for (int j=0;j<nxyz;j++) {float c0=freak.c[6*nxyz+j];if ((c0>4e-4)&&(c0<6e-4)) {std::cout << "cp "<<j<<" ";for (int i=0; i<ph_ncomp;i++){std::cout<< freak.c[i*nxyz+j]<<" ";} std::cout<<"\n";} }
	//for (int i=0;i<nxyz;i++) {std::cout<<"p "<<freak.p[i]<<" gvol "<<freak.gvol[i]<<"\n";}
	//writes to the dump file
	//bool dump_on = true;
	//bool append = false;
	//status = freak.PhreeqcRM_ptr->SetDumpFileName("phreeqc.dmp");
	//status = freak.PhreeqcRM_ptr->DumpModule(dump_on, append);
	/*if (ph_gcomp>0) {
	//std::cout<<"nb gcomp after "<<gcompn.size()<<" gvol 0 "<<freak.gvol[0]<<" sat "<<freak.wsat[0]<<"\n";
	std::cout<<"g 0 ";for (int i=0;i<ph_gcomp;i++) {std::cout <<freak.g[i*nxyz+0]<<" ";} std::cout<<"\n";
	std::cout<<"gm 0 ";for (int i=0;i<ph_gcomp;i++) {std::cout <<freak.gm[i*nxyz+0]<<" ";} std::cout<<"\n";
	std::cout<<"p ";for (int i=0;i<5;i++) {std::cout<<freak.p[i]<<" ";} std::cout<<"\n";
	}*/
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
	
int getSelOutput(my_phq & freak)
{
	try
	{
	IRM_RESULT status;
	//PhreeqcRM* phreeqc_rm = freak.PhreeqcRM_ptr;
	int n_user = freak.PhreeqcRM_ptr->GetNthSelectedOutputUserNumber(0);
	status = freak.PhreeqcRM_ptr->SetCurrentSelectedOutputUserNumber(n_user);
	int ncol = freak.PhreeqcRM_ptr->GetSelectedOutputColumnCount();std::cout<<"in phq ncol "<<ncol<<"\n";
	/*std::vector<std::string> headings(ncol);
	int nphase=0;
	for (int i = 0; i < ncol; i++)
		{
			status = freak.PhreeqcRM_ptr->GetSelectedOutputHeading(i,headings[i]);
			if (headings[i].substr(0,2)=="d_")  {nphase +=1;}
		}*/
	nxyz = freak.data[0];//std::cout<<"in phq nphase "<<nphase<<"\n";
	freak.nselect = ncol-6; //-nphase;
	freak.spc.resize((freak.nselect)*nxyz,0);
	std::vector<double> so(ncol*nxyz);
	status = freak.PhreeqcRM_ptr->GetSelectedOutput(so);
	for (int j = 0; j < nxyz; j++)
	{
		//nphase = 0; //reset nphase to 0 for reading below
		for (int i = 6; i < ncol; i++)
		{
			//if (headings[i].substr(0,2)=="d_") {nphase +=1;}
			//else 
				//{
			freak.spc[(i-6)*nxyz + j] = so[i*nxyz + j];
			if (j==10) {std::cout<< "selout sp "<< i << " : " << so[i*nxyz + j] << "\n";}
				//}
		}
	}
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


void register_basic_callback(void *cookie)
{
	my_phq *data;
#ifdef USE_MPI
	int mpi_tasks, mpi_myself;
#endif
	int	method_number = 1001;
	data = (my_phq *) cookie;

#ifdef USE_MPI
	MPI_Comm_size(MPI_COMM_WORLD, &mpi_tasks);
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_myself);
	if (mpi_myself == 0)
	{
		MPI_Bcast(&method_number, 1, MPI_INT, 0, MPI_COMM_WORLD);
	}
#endif

	const std::vector<IPhreeqcPhast *> w = data->PhreeqcRM_ptr->GetWorkers();
	for (int i = 0; i < (int) w.size(); i++)
	{
		w[i]->SetBasicCallback(my_basic_callback, cookie);
	}
}

double my_basic_callback(double x1, double x2, const char *str, void *cookie)
{
	my_phq * data_ptr = (my_phq *) cookie;
	PhreeqcRM *phreeqcrm_ptr = data_ptr->PhreeqcRM_ptr;
	std::string option(str);

	int cn = (int) x1;
	if (cn >= 0 && cn < phreeqcrm_ptr->GetChemistryCellCount())
	{
		//const std::vector < std::vector <int> > & back = phreeqcrm_ptr->GetBackwardMapping();
		if (option == "WSAT")
		{
			//int cn = back[rm_cell_number][0];std::cout<<cn<<" wsat "<<data_ptr->wsat[cn]<<"\n";
			return data_ptr->wsat[cn];
		}
		if (option == "PORO") { return data_ptr->poro[cn]; }
		if (option == "FOC") { return data_ptr->foc[cn]; }
	}
	return -999.9;
}
