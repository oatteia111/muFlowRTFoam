#include "PhreeqcRM.h" 
#include "IPhreeqc.hpp"
#include "IPhreeqcPhast.h"

#include <stdlib.h>
#include <vector>
#include <fstream>
#include <iterator>
#include <sstream>
#include <string>
#include <chrono> 

#include "partage.h"

double my_basic_callback(double x1, double x2, const char *str, void *cookie);
void register_basic_callback(void *cookie);

my_phq freak;

my_phq::my_phq()
{
	PhreeqcRM_ptr = NULL; //*PhreeqcRM_ptr;
#ifdef USE_MPI
	MPI_Comm rm_commxx;
#endif

}

int my_phq::phqInit() //my_phq & freak
{
	std::cout<<"in PHQINIT *********\n";
	nxyz = this->data[0];ph_ncomp = this->data[1];ph_gcomp = this->data[2];
	ph_nsolu=this->data[3];std::cout << " ph nxyz "<< nxyz <<"\n";
	int i,j;
	if (ph_ncomp<0) {bool EK = true;this->EK = EK;}
	else {std::vector<double> c(nxyz*ph_ncomp,0.);this->c = c;}
	
	#ifdef USE_MPI
		PhreeqcRM* phreeqc_rm = this->PhreeqcRM_ptr;
	#else
		
		int nthreads;
		PhreeqcRM *phreeqc_rm = new PhreeqcRM (nxyz, nthreads);
		this->PhreeqcRM_ptr = phreeqc_rm;		
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
		status = phreeqc_rm->SetUnitsSolution(2);           // 1 mg/L; 2 mol/L; 3 kg/kgs
		status = phreeqc_rm->SetUnitsPPassemblage(0);       // 0 mol/L cell; 1 mol/L water; 2 mol/L rock
		status = phreeqc_rm->SetUnitsExchange(0);           // 0 mol/L cell; 1 mol/L water; 2 mol/L rock
		status = phreeqc_rm->SetUnitsSurface(0);            // 0 mol/L RV; 1 mol/L water; 2 mol/L rock
		status = phreeqc_rm->SetUnitsGasPhase(0);           // 0 mol/L RV; 1 mol/L water; 2 mol/L rock
		status = phreeqc_rm->SetUnitsSSassemblage(1);       // 0 mol/L cell; 1 mol/L water; 2 mol/L rock
		status = phreeqc_rm->SetUnitsKinetics(0);           // 0 mol/L cell; 1 mol/L water; 2 mol/L rock
		if (this->EK) {status = phreeqc_rm->SetSpeciesSaveOn(true);}
				
		std::vector<double> rv(nxyz, 1.);
		status = phreeqc_rm->SetRepresentativeVolume(rv);   // RV is one dm3 of medium
		std::cout<<"poro size "<<this->poro.size()<<"\n";
		status = phreeqc_rm->SetPorosity(this->poro);
		this->deltaMvol.resize(nxyz,0.);

		// Open for writing file
		//std::quoted prf = "phq_out";
		status = phreeqc_rm->SetFilePrefix("debug");
		phreeqc_rm->OpenFiles();
		std::cout << "Outfile opened " << "\n";
		//set mask
		std::vector<int> print_chemistry_mask;
		print_chemistry_mask.resize(nxyz, 1);
		//for (int i=0;i<nxyz;i++) {print_chemistry_mask[i]=1;}
		status = phreeqc_rm->SetPrintChemistryMask(print_chemistry_mask);

		//load database and chemistry files and make the first run of the chem file
		status = phreeqc_rm->SetPrintChemistryOn(true, true, false); // workers, initial_phreeqc, utility
		std::cout << "dbase " << std::string(this->DB)<<"\n";
		status = phreeqc_rm->LoadDatabase(this->DB);
		std::cout << "dbase opened " << "\n";

		// Demonstrate add to Basic: Set a function for Basic CALLBACK after LoadDatabase
		register_basic_callback(&freak);

		bool workers = true;             // Worker instances do the reaction calculations for transport
		bool initial_phreeqc = true;     // InitialPhreeqc instance accumulates initial and boundary conditions
		bool utility = true;             // Utility instance is available for processing
		status = phreeqc_rm->RunFile(workers, initial_phreeqc,utility, this->ChemFile);

		// Clear contents of workers and utility
		initial_phreeqc = false;
		std::string input = "DELETE; -all";
		status = phreeqc_rm->RunString(workers, initial_phreeqc, utility, input.c_str());
		
		int ph_ncomp = phreeqc_rm->FindComponents(); // NECESSARY don't remove
		const std::vector<std::string> &gcomp = phreeqc_rm->GetGasComponents(); //std::cout<<"phq gcomp0 "<<gcomp[0]<<"\n";
		const std::vector<std::string> &comp = phreeqc_rm->GetComponents(); 
		this->comp = comp;this->gcomp = gcomp;
		std::cout<<"phq comp "; for (i=0;i<ph_ncomp;i++) {this->comp[i] = comp[i];std::cout<<comp[i]<<" ";} std::cout <<"\n";

		// Set array of initial conditions by reading the solution... numbers in data
		std::vector<int> ic1;
		ic1.resize(nxyz*7, -1);std::cout<<"ic1 size "<<ic1.size()<<"\n";
		//std::vector<double> f1;
		int icnt = 4;int mx; // to consider the first 5 numbers : nxyz,ncomp,gcomp,nsolu,units (++ in the first line of the for makes 5 numbers)
		for (i = 0; i < 7; i++) // Solution,Equilibrium phases,Exchange,Surface,Gas phase,Solid solutions,Kinetics
		{
			icnt++;mx=0;
			int test = this->data[icnt];//std::cout<<"index "<<i<<" " <<test<<"\n";
			if (test == 0)  // if 0 just one value to read
				{
				icnt++; 
				for (int j=0;j<nxyz;j++) {
					ic1[i*nxyz+j] = this->data[icnt];
					//if ((j==0)||(j==nxyz-2)) {std::cout << i << " cell "<< j << " " << this->data[icnt]<< "\n";} 
					}
				}
			else 		
				{
				for (int j=0;j<nxyz;j++) 
					{
					icnt++;ic1[i*nxyz+j] = this->data[icnt]; 
					mx = std::max(mx,this->data[icnt]);
					//if ((j==0)||(j==nxyz-2))  {std::cout << i << " cell "<< j << " " << this->data[icnt]<< "\n";} 
					}
				}
			if (i==1) {this->nphase=mx;}
			if (i==2) {this->nexch=mx;}
			if (i==3) {this->nsurf=mx;}
		}
		std::cout << "end fill here\n";
		/*
		const std::vector<IPhreeqcPhast *> w = phreeqc_rm->GetWorkers();
		w[0]->AccumulateLine("Delete; -all");
		int w0 = w[0]->RunAccumulated();
		*/
		//set initial conditions
		status = phreeqc_rm->InitialPhreeqc2Module(ic1);
		std::cout << "end init Ic \n";
		
		// run initial conditions
		int ncell = phreeqc_rm->GetChemistryCellCount();
		std::cout << "ncell "<< ncell << " ph_ncomp "<< ph_ncomp<< " starting init phreeqc \n";
		status = phreeqc_rm->SetTemperature(this->temp);std::cout<<"temp "<<this->temp[0]<<"\n";
		if (ph_gcomp>0) {
			//status = this->PhreeqcRM_ptr->SetSaturation(this->wsat); // keep the order of three lines to set values
			//status = phreeqc_rm->SetPressure(this->p);
			status = this->PhreeqcRM_ptr->SetGasPhaseVolume(this->gvol);
			status = this->PhreeqcRM_ptr->SetGasCompMoles(this->gm);
		//std::cout<<"gvol "; for (j=0;j<5;j++){std::cout<<this->gvol[j]<<" ";}std::cout<<"\n";
		//std::cout<<"gm "; for (j=0;j<5;j++){ for (i=0;i<ph_gcomp;i++) {std::cout<<this->gm[i*nxyz+j]<<" ";} std::cout<<"\n";}
			}
		status = phreeqc_rm->RunCells();
		
		std::cout << " cells run \n";
		status = phreeqc_rm->GetConcentrations(this->c);
		
		//set boundary condition
		std::vector<int> bc1;
		bc1.resize(ph_nsolu,0);
		for (i=0;i<ph_nsolu;i++) {bc1[i]=i;}//std::cout<<"in phq bc "<<i<<"\n";}
		status = phreeqc_rm->InitialPhreeqc2Concentrations(this->bc_conc, bc1);
		std::cout<<" bc "<<this->bc_conc[0]<<"\n";
		if (this->EK>0) {
			int nspc = phreeqc_rm->GetSpeciesCount(); std::cout<<"phq eK nspc "<<nspc<<"\n";
			this->nspc = nspc;
			const std::vector<std::string> &spcNames = phreeqc_rm->GetSpeciesNames();
			//for (i=0;i<this->nspc;i++) {std::cout<<" "<<spcNames[i];} std::cout<<"\n";
			const std::vector <double> &z = phreeqc_rm->GetSpeciesZ();
			this->z = z;
			//for (i=0;i<nspc;i++) {std::cout<<z[i]<<"\n";}
			std::vector<double> c(nxyz*nspc,0.);this->c = c;
			status = phreeqc_rm->GetSpeciesConcentrations(this->c);
			status = phreeqc_rm->InitialPhreeqc2SpeciesConcentrations(this->bc_species, bc1);
			//std::cout<<"phq spc \n"; 
			//for (j=0;j<nxyz;j++){for (i=0;i<nspc;i++) {if (this->z[i]!=0){std::cout<<this->c[i*nxyz+j]<<" ";}} std::cout<<"\n";}
			const std::vector <double> &diff25 = phreeqc_rm->GetSpeciesD25();
			this->diff25 = diff25; 
			//std::cout<<"phq eK diff25 0 ";for (i=0;i<nspc;i++) {std::cout<<" "<<diff25[i]<<"\n";}
		}
		if (ph_gcomp>0) {
			status = phreeqc_rm->GetSaturation(this->wsat);
			status = phreeqc_rm->GetGasCompPressures(this->g);
			status = phreeqc_rm->GetGasCompMoles(this->gm);
			//status = phreeqc_rm->GetGasPhaseVolume(this->gvol);
			//std::cout << "in phq poro "<<this->poro[1]<<"\n";
			//std::cout << "c 1-4 ";for (int i=0; i<4;i++){std::cout<< this->c[i*nxyz+1]<<" ";} std::cout<<"\n";
			//std::cout << "c 5-9 ";for (int i=4; i<8;i++){std::cout<< this->c[i*nxyz+1]<<" ";} std::cout<<"\n";
			//std::cout<<"Vm "<<this->gvol[0]/this->gm[0]<<" "<<this->gvol[0]/this->gm[1]<<"\n"; // false
		//std::cout<<"gvol "; for (j=0;j<5;j++){std::cout<<this->gvol[j]<<" ";}std::cout<<"\n";
		//std::cout<<"gm "; for (j=0;j<5;j++){ for (i=0;i<ph_gcomp;i++) {std::cout<<this->gm[i*nxyz+j]<<" ";} std::cout<<"\n";}
			for (i=0;i<gcomp.size();i++) { if (gcomp[i]=="H2O(g)") {this->iGwater = i;} }
			//std::cout << "c 0 ";for (int i=0; i<ph_ncomp;i++){std::cout<< this->c[i*nxyz+0]<<" ";} std::cout<<"\n";
		}
	//std::cout << "c 1 ";for (int i=0; i<ph_ncomp;i++){std::cout<< this->c[i*nxyz+1]<<" ";} std::cout<<"\n";
	std::cout <<"end phq init "<<"\n";
		
		//printing to file
		status = phreeqc_rm->CloseFiles();
		//status = this->PhreeqcRM_ptr->SetDumpFileName("phqAll.dmp");
		//status = phreeqc_rm->SetPrintChemistryOn(true, false, true); // workers, initial_phreeqc, utility
		// below : trial, does not work	
		//const std::vector<std::string> &PPcomp = this->PhreeqcRM_ptr->phreeqc_bin->Get_PPassemblage(0)->cxxPPassemblage->Get_EltList();
		//std::cout<<"PP comp "; for (i=0;i<PPcomp.size();i++) {std::cout<<PPcomp[i]<<" ";} std::cout <<"\n";

	
	}
	catch (PhreeqcRMStop)
	{
		std::cerr << "failed with an error in PhreeqcRM." << std::endl;
	/* we don't want to have mpi abort as we can now handle error son phreeqc
	#ifdef USE_MPI
		MPI_Abort(MPI_COMM_WORLD, 1);
	#endif
	*/
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

int my_phq::phqRun() //my_phq & freak
{
	try
	{
	IRM_RESULT status;
	//PhreeqcRM* phreeqc_rm = this->PhreeqcRM_ptr;
	status = this->PhreeqcRM_ptr->SetErrorHandlerMode(1);
	//status = this->PhreeqcRM_ptr->SetFilePrefix("phq_run");
	//this->PhreeqcRM_ptr->OpenFiles();
	//std::vector<int> print_chemistry_mask;
	//print_chemistry_mask.resize(nxyz, 1);//print_chemistry_mask[0]=1;print_chemistry_mask[50]=1;
	//status = this->PhreeqcRM_ptr->SetPrintChemistryMask(print_chemistry_mask); 
	//status = this->PhreeqcRM_ptr->SetPrintChemistryOn(true, false, false); // workers, initial_phreeqc, utility
	const std::vector<std::string> &components = this->PhreeqcRM_ptr->GetComponents();
	const std::vector<std::string> &gcomp = this->PhreeqcRM_ptr->GetGasComponents(); std::cout<<"phq ngcomp sz "<<gcomp.size()<<" phq ws "<<this->wsat[0]<<"\n";
	nxyz = this->data[0];
	std::cout << "in phq nxyz "<<nxyz<<" size c "<<this->c.size()<<" poro "<<this->poro[2]<<"\n";
	
	status = this->PhreeqcRM_ptr->SetTemperature(this->temp);
	//std::cout<<"temp 0-5 ";for (int j=0; j<5;j++){std::cout<<this->temp[j]<<" ";} std::cout<<"\n";
	//std::cout << "c 0 ";
	//for (int j=0;j<8;j++) {for (int i=0; i<ph_ncomp;i++){std::cout<< this->c[i*nxyz]<<" ";} std::cout<<"\n";}
	//std::cout << "c mid ";for (int i=0; i<ph_ncomp;i++){std::cout<< this->c[i*nxyz+int(nxyz/2)]<<" ";} std::cout<<"\n";
	//std::cout << "c end ";for (int i=0; i<ph_ncomp;i++){std::cout<< this->c[i*nxyz+nxyz-2]<<" ";} std::cout<<"\n";
	//std::cout<<"press 0-5 ";for (int j=0; j<5;j++){std::cout<<this->p[j]<<" ";} std::cout<<"\n";//don't use it pressure is calculated by phq
	//for (int j=0;j<nxyz;j++) {float c0=this->c[6*nxyz+j];if ((c0>4e-4)&&(c0<6e-4)) {std::cout << "cp "<<j<<" ";for (int i=0; i<ph_ncomp;i++){std::cout<< this->c[i*nxyz+j]<<" ";} std::cout<<"\n";} }
	if (gcomp.size()>0) {
		//std::cout<<"nb gcomp before "<<gcomp.size()<<"\n";//<<" gvol 3 "<<this->gvol[0]<<" sat "<<this->wsat[0]<<"\n";
		//std::cout<<"g mid ";for (int i=0;i<ph_gcomp;i++) {std::cout <<this->g[i*nxyz+int(nxyz/2)]<<" ";} std::cout<<"\n";
		//std::cout<<"gm mid ";for (int i=0;i<ph_gcomp;i++) {std::cout <<this->gm[i*nxyz+int(nxyz/2)]<<" ";} std::cout<<"\n";
		status = this->PhreeqcRM_ptr->SetSaturation(this->wsat); // keep the 3 lines and the order of three lines to set values
		status = this->PhreeqcRM_ptr->SetGasPhaseVolume(this->gvol);
		status = this->PhreeqcRM_ptr->SetGasCompMoles(this->gm);
		//status = this->PhreeqcRM_ptr->SetPressure(this->p); // we cannot set preesur ein fact
		//for (int j=0;j<8;j++) {std::cout<<j<<" p "<<this->p[j]<<" gvol "<<this->gvol[j]<<" wsat "<<this->wsat[j]<<" g1 "<<this->gm[nxyz+j]<<"\n";}
		}
	else 
		{
		status = this->PhreeqcRM_ptr->SetSaturation(this->wsat);
		}
	//std::cout<<"phq wsat[-2] "<<this->wsat[nxyz-2];
	status = this->PhreeqcRM_ptr->GetSaturation(this->wsat);
	//std::cout<<"phq wsat[-2] "<<this->wsat[nxyz-2];

	if (this->EK) {status = this->PhreeqcRM_ptr->SpeciesConcentrations2Module(this->c);}
	else {status = this->PhreeqcRM_ptr->SetConcentrations(this->c);}
	status = this->PhreeqcRM_ptr->SetTimeStep(this->tstep);std::cout<<"tstep "<<this->tstep<<"\n";
	//----- run here ------
	status = this->PhreeqcRM_ptr->RunCells();
	std::cout<<"run done \n";
	
	//status = this->PhreeqcRM_ptr->GetSaturation(this->wsat);
	//status = this->PhreeqcRM_ptr->GetTemperature(this->temp);
	if (this->EK) {status = this->PhreeqcRM_ptr->GetSpeciesConcentrations(this->c);}
	else {status = this->PhreeqcRM_ptr->GetConcentrations(this->c);}
	
	//for (int j=0;j<8;j++) {for (int i=0; i<ph_ncomp;i++){std::cout<< this->c[i*nxyz]<<" ";} std::cout<<"\n";}
	
	if (gcomp.size()>0) {
		status = this->PhreeqcRM_ptr->GetGasCompPressures(this->g);
		status = this->PhreeqcRM_ptr->GetGasCompMoles(this->gm);
		status = this->PhreeqcRM_ptr->GetGasPhaseVolume(this->gvol);
		this->p = this->PhreeqcRM_ptr->GetPressure();
		//for (int j=0;j<8;j++) {std::cout<<j<<" p "<<this->p[j]<<" gvol "<<this->gvol[j]<<" ";for (int i=0;i<gcomp.size();i++) {std::cout<<this->g[i*nxyz+j]<<" "<<this->gm[i*nxyz+j]<<" ";}std::cout<<"\n";}
		//std::cout<<"nb gcomp after "<<ph_gcomp<<" gvol mid "<<this->gvol[int(nxyz/2)]<<" sat "<<this->wsat[int(nxyz/2)]<<"\n";
		/*
		std::cout<<"g mid ";for (int i=0;i<ph_gcomp;i++) {std::cout <<this->g[i*nxyz+int(nxyz/2)]<<" ";} std::cout<<"\n";
		std::cout<<"gm mid ";for (int i=0;i<ph_gcomp;i++) {std::cout <<this->gm[i*nxyz+int(nxyz/2)]<<" ";} std::cout<<"\n";
		std::cout<<"p mid ";std::cout<<this->p[nxyz/2]<<"\n";
		*/
		}
	std::cout<<"end phq run \n ";
	}
	catch (PhreeqcRMStop)
	{
		std::cerr << "failed with an error in PhreeqcRM." << std::endl;
	/*
	#ifdef USE_MPI
		MPI_Abort(MPI_COMM_WORLD, 1);
	#endif
	*/
		return IRM_FAIL;
	}
	catch (...)
	{
		std::cerr << "failed with an unhandled exception." << std::endl;
	/*
	#ifdef USE_MPI
		MPI_Abort(MPI_COMM_WORLD, 1);
	#endif
	*/
		return IRM_FAIL;
	}
	return EXIT_SUCCESS;
};
	
int my_phq::getSelOutput() //my_phq & freak
{
	try
	{
	IRM_RESULT status;
	//PhreeqcRM* phreeqc_rm = this->PhreeqcRM_ptr;
	//int n_user = this->PhreeqcRM_ptr->GetNthSelectedOutputUserNumber(0);
	//status = this->PhreeqcRM_ptr->SetCurrentSelectedOutputUserNumber(n_user);
	int ncol = this->PhreeqcRM_ptr->GetSelectedOutputColumnCount();//std::cout<<"in phq ncol "<<ncol<<"\n";
	std::vector<std::string> headings(ncol);
	//int nphase=0;
	for (int i = 0; i < ncol; i++)
		{
			status = this->PhreeqcRM_ptr->GetSelectedOutputHeading(i,headings[i]);
			//if (headings[i].substr(0,2)=="d_")  {nphase +=1;}
		}
	nxyz = this->data[0];//std::cout<<"in phq selout nxyz "<<nxyz<<" for cell 0 and nxyz-2\n";
	std::vector<double> so(ncol*nxyz);
	status = this->PhreeqcRM_ptr->GetSelectedOutput(so);
	this->nselect = ncol-6; //-nphase;
	this->spc.resize((this->nselect)*nxyz,0);
	for (int j = 0; j < nxyz; j++)
	{
		if (this->wsat[j]>0) {
			for (int i = 6; i < ncol; i++)
			{
				if (headings[i].substr(0,2)=="d_") {this->deltaMvol[j]+=so[i*nxyz + j];}
				this->spc[(i-6)*nxyz + j] = so[i*nxyz + j];
				//if ((j==0)||(j==nxyz-2)) {std::cout<< "selout sp "<< i << " : " << so[i*nxyz + j] << "\n";}
					//}
			}
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

int my_phq::phqTest(std::string time) {
	//IRM_RESULT status;
	//int ph_ncomp = this->PhreeqcRM_ptr->FindComponents();
	return 1;
	};

int my_phq::phqRestart(std::vector<int> & ph_data, std::string time,std::vector<double> c_ph,std::vector<double> g_ph) 
{
	//namespace fs = std::filesystem;
	IRM_RESULT status;int ncomp,ngcomp,nphas;std::cout<<"in restart \n";
	//read data in the given timestep t
	std::ifstream inputSpc{cur_dir+"/time/Species"};
	//std::ofstream spc2{cur_dir/time/"Species_bak"}; 
	//std::string line;
	//while (std::getline(inputSpc, line)) {spc2 << line << "\n"; } 
	std::vector<double>Species = {std::istream_iterator<float>{inputSpc}, {}};

	// the idea is to build a long string with the instructions and to runString (but phreeqc ha not been run before)

	//nxyz = this->data[0];
	//int ph_ncomp = this->PhreeqcRM_ptr->FindComponents(); 
	const std::vector<std::string> comp = this->PhreeqcRM_ptr->GetComponents();// std::cout<<"phq comp0 "<<comp[0]<<"\n";
	const std::vector<std::string> gcomp = this->PhreeqcRM_ptr->GetGasComponents();																				
	const std::vector<std::string> phasStr = this->PhreeqcRM_ptr->GetEquilibriumPhases(); 
	int nspec = this->PhreeqcRM_ptr->GetSelectedOutputColumnCount()-6;
	nphas = phasStr.size();
	//get phqfoam data
	nxyz=ph_data[0];ncomp=ph_data[1];ngcomp=ph_data[2]; //!!! nxyz here is inside the ractive part
	//read old initchem to get the list of exchange & surface (one exch or surf assemblage/line)
	std::vector<std::string> initChem;
	std::ifstream inputChem{cur_dir+"/initChem.pqi" };
	initChem = {std::istream_iterator<std::string>{inputChem}, {}};
	std::vector<std::string> seloStr,exchStr,surfStr,gasStr,kinStr;
	int flagSelo,flagE,flagS,flagG,flagK;
	flagSelo=1;int cnt,exlength,sulength,glength;
	for (size_t i=0;i<initChem.size();i++)
		{
		cnt+=1;
		if (flagSelo==1) {seloStr.push_back(initChem[i]);}
		if (initChem[i]=="Solution") {flagSelo=0;}
		if (initChem[i]=="Exchange") {flagE=1;cnt=0;}
		if (initChem[i]=="Surface") {flagS=1;cnt=0;}
		if (initChem[i]=="-equilibrate") {
			if (flagE==1) {exlength=cnt;}
			if (flagS==1) {sulength=cnt;}		   
			flagE=0;flagS=0;flagG=0;cnt=0;
			}
		if (initChem[i]=="Kinetics") //supposed to be at the end of the file
			{flagK=1;cnt=0;
			for (int j=i;j<initChem.size();j++){kinStr.push_back(initChem[j]);}
			break;
			}
		if (flagE==1) {exchStr.push_back(initChem[i]);}
		if (flagS==1) {surfStr.push_back(initChem[i]);}
		}
	int nexch=this->nexch;	int nsurf=this->nsurf;
	std::cout<<"inichem read , nexch, nsurf "<< nexch<<" "<<nsurf<<" n data "<<ph_data.size()<<"\n";
	//Solution,Equilibrium phases,Exchange,Surface,Gas phase,Solid solutions,Kinetics
	std::ofstream outPhq{cur_dir+"/phqfoam1.txt"};
	for (int i=0;i<5;i++) {outPhq<<ph_data[i]<<" ";} outPhq<<"\n";
	outPhq<<"-1 ";for (int j=0;j<nxyz;j++) {outPhq<<j<<" ";} outPhq<<"\n";  //solutions
	if (nphas>0) {outPhq<<"-1 ";for (int j=0;j<nxyz;j++) {outPhq<<j<<" ";} outPhq<<"\n";} //phases
	else {outPhq<<"0 -1 \n";}
	if (nexch>0) {outPhq<<"-1 ";for (int j=0;j<nxyz;j++) {outPhq<<j<<" ";} outPhq<<"\n";}
	else {outPhq<<"0 -1 \n";}
	if (nsurf>0) {outPhq<<"-1 ";for (int j=0;j<nxyz;j++) {outPhq<<j<<" ";} outPhq<<"\n";}
	else {outPhq<<"0 -1 \n";}	
	if (ngcomp>0) {outPhq<<"-1 ";for (int j=0;j<nxyz;j++) {outPhq<<j<<" ";} outPhq<<"\n";}
	else {outPhq<<"0 -1 \n";}						   
	outPhq<<"0 -1 \n";
	if (flagK==1) {outPhq<<"0 1 \n";} else {outPhq<<"0 -1 \n";}
	outPhq.close();
		
	//writes the new initchem
	int ncol = this->PhreeqcRM_ptr->GetSelectedOutputColumnCount();
	std::vector<std::string> headings(ncol);
	int k=0;std::vector<int> indxP(nphas);
	if (nphas>0)	
	{
		for (int i = 0; i < ncol; i++)
		{
			status = this->PhreeqcRM_ptr->GetSelectedOutputHeading(i,headings[i]);
			if (headings[i]==phasStr[k])  {
				indxP[k]=i-6;std::cout<<" phas "<<i<<" "<<k<<" "<<headings[i]<<" "<<phasStr[k]<<" "<<indxP[k]<<"\n";
				k+=1;
			}
		}
	}
	std::ofstream outChem{cur_dir+"/initChem1.pqi"};
	//std::string oo = "-";
	int nexchj,nsurfj,ngasj,start;
	for (size_t i=0;i<seloStr.size()-1;i++) { if (seloStr[i][0]=='-') {outChem<<";";} outChem<<seloStr[i]<<" ";} outChem<<"\n";
	for (int j=0;j<nxyz;j++)
		{
		outChem<<"Solution "<<j<<"; units mol/kgw ;"; // find the solution (this has been read by opf) st
		for (int i=4;i<ncomp;i++) {outChem << comp[i]<<" "<<c_ph[i*nxyz+j]<<";"; } // conc in solution
		//outChem <<" pH "<<Species[j*nspec]<<"; pe "<<Species[j*nspec+1]<<"; ";
		if (nphas>0) {outChem<<"Equilibrium_Phases "<<j<<"; ";}
		for (int i=0;i<nphas;i++) {outChem<<phasStr[i]<<" 0 "<<Species[j*nspec+indxP[i]]<<";";}
		if (nexch>0) {
			outChem<<"Exchange "<<j<<"; ";nexchj=ph_data[4+(nxyz+1)*2+1+j];start=nexchj*(exlength-2);
			for (int i=start;i<start+exlength-2;i++) {outChem<<exchStr[i]<<" ";}
			outChem<<";-e "<<j<<";";
			}
		if (nsurf>0) {
			outChem<<"Surface "<<j<<"; ";nsurfj=ph_data[4+(nxyz+1)*3+1+j];start=nsurfj*sulength;
			for (int i=start;i<start+sulength;i++) {outChem<<surfStr[i]<<" ";}
			outChem<<";-e "<<j<<";-no_edl ";
			}
		if (ngcomp>0) {
			outChem<<"Gas_Phase "<<j<<";-fixed_pressure;-pressure 1;";// find the gas (this has been read by opf) st
			for (int i=0;i<ngcomp;i++) {outChem << gcomp[i]<<" "<<g_ph[i*nxyz+j]/this->gvol[j]<<";"; } // conc in solution
			}
		outChem<<"\n";
		}  
	if (flagK==1) 
		{
		for (int i=0;i<kinStr.size();i++) {outChem<<kinStr[i]<<" ";}
		}
	outChem.close();	
	/*
	std::string inputPhq;//std::cout<<" test "<<c_ph[4*nxyz+233]<<"\n";
	for (int j=0;j<nxyz;j++)
		{
		string js=std::to_string(j);
		inputPhq+="Solution "+js+"; units mol/kgw ;"; // find the solution (this has been read by opf) st
		for (int i=4;i<ncomp;i++) {inputPhq += comp[i]+" "+std::to_string(c_ph[i*nxyz+j])+";"; } // conc in solution
		//inputPhq <<" pH "<<Species[j*nspec]<<"; pe "<<Species[j*nspec+1]<<"; ";
		inputPhq+="Equilibrium_Phases "+js+"; ";
		for (int i=0;i<nphas;i++) {inputPhq+=phasStr[i]+" 0 "+std::to_string(Species[j*nspec+2+i*2])+";";}
		inputPhq+="Exchange "+js+"; "; // just write equilibrate with exchanger
		for (int i=2;i<nexch;i++) {inputPhq+=exchStr[i]+" ";}
		inputPhq+=";-e "+js+";";
		inputPhq+="Surface "+js+"; "; // just write equilibrate with surface
		for (int i=2;i<nsurf;i++) {inputPhq+=surfStr[i]+" ";}
		inputPhq+=";-e "+js+";-no_edl \n";
		}  
	std::cout<<"end string \n";
	std::string inS="Equilibrium_Phases 10;Fluorite 0 1";
	status = this->PhreeqcRM_ptr->RunString(false, true, false, inS.c_str());  // bool : workers, initial, utility
	status = this->PhreeqcRM_ptr->RunString(false, true, false, inputPhq.c_str());  // bool : workers, initial, utility
	status = this->PhreeqcRM_ptr->RunFile(false, true, false, cur_dir+"/restart.txt");  // bool : workers, initial, utility
	std::cout<<"phq string run \n";
	*/
	return 1;
};


void my_phq::register_basic_callback(void *cookie)
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
		w[i]->SetBasicCallback(&my_phq::my_basic_callback, cookie);
	}
}

double my_phq::my_basic_callback(double x1, double x2, const char *str, void *cookie)
{
	my_phq *data_ptr = (my_phq *) cookie;
	PhreeqcRM *phreeqcrm_ptr = data_ptr->PhreeqcRM_ptr;
	std::string option(str);

	int cn = (int) x1;
	if (cn >= 0 && cn < phreeqcrm_ptr->GetChemistryCellCount())
	{
		const std::vector < std::vector <int> > & back = phreeqcrm_ptr->GetBackwardMapping();
		if (option == "WSAT")
		{
			//int cn = back[rm_cell_number][0];
			//std::cout<<cn<<" wsat "<<data_ptr->wsat[cn]<<"\n";
			return data_ptr->wsat[cn];
		}
		if (option == "PORO") { return data_ptr->poro[cn]; }
		if (option == "FOC")  { return data_ptr->foc[cn]; }
		if (option == "TEMP") { return data_ptr->temp[cn]; }
	}
	return -999.9;
}
