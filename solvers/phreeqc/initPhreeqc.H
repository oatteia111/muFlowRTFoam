#if defined(USE_MPI)
#include <mpi.h>
#endif
#include "PhreeqcRM.h"
#include "IPhreeqc.hpp"
#include "IPhreeqcPhast.h"
int nthreads = 3;

///////////////// get some important values
/* ifstream f0(cur_dir/"phqfoam.txt");
std::vector<int> ph_data;
if (f0) {        
    int value;
    // read the elements in the file into a vector  
    while ( f0 >> value ) {
        ph_data.push_back(value);
    }
f0.close(); */

int nxyz,ph_ncomp,ph_gcomp,ph_nsolu,sunits;
////////////////////// start phreeqc
class my_phq
{
public:
	PhreeqcRM *PhreeqcRM_ptr;
	int nxyz=3;
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

	//std::vector<std::string> complist;
    //void setCompList(std::vector<std::string>) {this->complist = complist;}
	//this->complist = complist;
	void init()
	{
		nxyz = this->data[0];ph_ncomp = this->data[1];ph_gcomp = this->data[2];sunits = this->data[4];std::cerr << " ph nxyz "<< nxyz <<"\n";
		int i,j;
		if (ph_ncomp<0) {bool EK = true;this->EK = EK;}
		else {std::vector<double> c(nxyz*ph_ncomp,0.);this->c = c;}
	/*#ifdef USE_MPI
		// MPI
		PhreeqcRM *phreeqc_rm = new PhreeqcRM(nxyz, MPI_COMM_WORLD);
		MP_TYPE comm = MPI_COMM_WORLD;
		int mpi_myself;
		if (MPI_Comm_rank(MPI_COMM_WORLD, &mpi_myself) != MPI_SUCCESS)
		{
			exit(4);
		}
		if (mpi_myself > 0)
		{
			phreeqc_rm.SetMpiWorkerCallbackC(worker_tasks_cc);
			phreeqc_rm.SetMpiWorkerCallbackCookie(&some_data);
			phreeqc_rm.MpiWorker();
			return EXIT_SUCCESS;
		}
	#else
		PhreeqcRM *phreeqc_rm = new PhreeqcRM (nxyz, nthreads);
	#endif */
		PhreeqcRM *phreeqc_rm = new PhreeqcRM (nxyz, nthreads);
		this->PhreeqcRM_ptr = phreeqc_rm;
		IRM_RESULT status;
		status = phreeqc_rm->SetErrorHandlerMode(1);
		status = phreeqc_rm->SetScreenOn(false);
		status = phreeqc_rm->SetComponentH2O(true); // to get water,excess H an dexcess O transported
		status = phreeqc_rm->SetRebalanceFraction(0.5);
		status = phreeqc_rm->SetRebalanceByCell(true);
		phreeqc_rm->UseSolutionDensityVolume(false);
		phreeqc_rm->SetPartitionUZSolids(false);
		status = phreeqc_rm->SetUnitsSolution(2);           // 1, mg/L; 2, mol/L; 3, kg/kgs
		status = phreeqc_rm->SetUnitsPPassemblage(sunits);       // 0, mol/L cell; 1, mol/L water; 2 mol/L rock
		status = phreeqc_rm->SetUnitsExchange(sunits);           // 0, mol/L cell; 1, mol/L water; 2 mol/L rock
		status = phreeqc_rm->SetUnitsSurface(sunits);            // 0, mol/L RV; 1, mol/L water; 2 mol/L rock
		status = phreeqc_rm->SetUnitsGasPhase(0);           // 0, mol/L RV; 1, mol/L water; 2 mol/L rock
		status = phreeqc_rm->SetUnitsSSassemblage(sunits);       // 0, mol/L cell; 1, mol/L water; 2 mol/L rock
		status = phreeqc_rm->SetUnitsKinetics(0);           // 0, mol/L cell; 1, mol/L water; 2 mol/L rock
		if (this->EK) {status = phreeqc_rm->SetSpeciesSaveOn(true);}
				
		std::vector<double> rv(nxyz, 1.);
		status = phreeqc_rm->SetRepresentativeVolume(rv);   // RV is one dm3 of medium
		std::cout<<"poro size "<<this->poro.size()<<"\n";
		status = phreeqc_rm->SetPorosity(this->poro);

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
		status = phreeqc_rm->LoadDatabase(this->DB);
		std::cout << "dbase opened " << "\n";
		bool workers = true;             // Worker instances do the reaction calculations for transport
		bool initial_phreeqc = true;     // InitialPhreeqc instance accumulates initial and boundary conditions
		bool utility = true;             // Utility instance is available for processing
		status = phreeqc_rm->RunFile(workers, initial_phreeqc,utility, this->ChemFile);

		// Clear contents of workers and utility
		initial_phreeqc = false;
		std::string input = "DELETE; -all";
		status = phreeqc_rm->RunString(workers, initial_phreeqc, utility, input.c_str());
		
		int ph_ncomp = phreeqc_rm->FindComponents(); // NECESSARY don't remove
		const std::vector<std::string> &gcompn = phreeqc_rm->GetGasComponents(); //std::cout<<"phq gcomp0 "<<gcomp[0]<<"\n";
		const std::vector<std::string> &compn = phreeqc_rm->GetComponents(); 
		std::cout<<"phq comp "; for (i=0;i<ph_ncomp;i++) {std::cout<<compn[i]<<" ";} std::cout <<"\n";

		// Set array of initial conditions by reading the solution... numbers in data
		std::vector<int> ic1, ic2;
		ic1.resize(nxyz*7, -1);
		//std::vector<double> f1;
		int icnt = 4; // to consider the first 5 numbers : nxyz,ncomp,gcomp,nsolu,units (++ in the first line of the for makes 5 numbers)
		for (i = 0; i < 7; i++) // Solution,Equilibrium phases,Exchange,Surface,Gas phase,Solid solutions,Kinetics
		{
			icnt++;
			int test = this->data[icnt];//std::cout<<"index "<<i<<" " <<test<<"\n";
			if (test == 0)  // if 0 just one value to read
				{icnt++; 
				for (int j=0;j<nxyz;j++) {ic1[i*nxyz+j] = this->data[icnt];
				if (j==10) {std::cout << i << " cell "<< j << " " << this->data[icnt]<< "\n";} }
				}
			else 		
				{for (int j=0;j<nxyz;j++) {icnt++;ic1[i*nxyz+j] = this->data[icnt];
				if (j==10) {std::cout << i << " cell "<< j << " " << this->data[icnt]<< "\n";} }
			}
		}
		const std::vector<IPhreeqcPhast *> w = phreeqc_rm->GetWorkers();
		w[0]->AccumulateLine("Delete; -all");
		int w0 = w[0]->RunAccumulated();

		status = phreeqc_rm->InitialPhreeqc2Module(ic1);
		
		// Set solutions composition at fixed cells HERE only TWO solutions (shall be modified)
		//std::vector<double> solu_conc;
		//std::vector<int> cf1 ={0,-1};
		//for (i= 0;i<2;i++) {cf1[this->cfix[2+i*4+2]] = this->cfix[2+i*4+3];std::cout<<" cfix "<<i<<" "<< this->cfix[2+i*4+2] << "\n";}

		// run initial conditions
		int ncell = phreeqc_rm->GetChemistryCellCount();
		std::cout << "ncell "<< ncell << " ph_ncomp "<< ph_ncomp<< " starting init phreeqc \n";
		//status = phreeqc_rm->SetSaturation(this->wsat);
		//status = phreeqc_rm->SetGasPhaseVolume(this->gvol);
		//status = phreeqc_rm->SetGasCompMoles(this->gm);
		//status = phreeqc_rm->SetPressure(this->p);
		//phreeqc_rm->SetPrintChemistryOn(true,false,false);
		status = phreeqc_rm->SetTemperature(this->temp);
		
		status = phreeqc_rm->RunCells();
		
		std::cout << " cells run \n";
		status = phreeqc_rm->GetConcentrations(this->c);
		if (this->EK) {
			int nspc = phreeqc_rm->GetSpeciesCount(); std::cout<<"phq eK nspc "<<EK<<" "<<nspc<<"\n";
			this->nspc = nspc;
			status = phreeqc_rm->GetSpeciesConcentrations(this->spc);
			const std::vector <double> & diff25 = phreeqc_rm->GetSpeciesD25();
			this->diff25 = diff25;
			const std::vector <double> & z = phreeqc_rm->GetSpeciesZ();
			this->z = z;
			const std::vector<std::string> & spcNames = phreeqc_rm->GetSpeciesNames();
			this->spcNames = spcNames;
		}
		if (ph_gcomp>0) {
			status = phreeqc_rm->GetSaturation(this->wsat);
			status = phreeqc_rm->GetGasCompPressures(this->g);
			status = phreeqc_rm->GetGasCompMoles(this->gm);
			status = phreeqc_rm->GetGasPhaseVolume(this->gvol);
			std::cout << "in phq poro "<<this->poro[1]<<"\n";
			std::cout << "c 1-4 ";for (int i=0; i<4;i++){std::cout<< this->c[i*nxyz+1]<<" ";} std::cout<<"\n";
			std::cout << "c 5-9 ";for (int i=4; i<8;i++){std::cout<< this->c[i*nxyz+1]<<" ";} std::cout<<"\n";
			//std::cout<<"gvol 0 "<<this->gvol[0]<<"gvol 1 "<<this->gvol[1]<<"\n";
			//std::cout<<"g 0 ";for (int i=0;i<ph_gcomp;i++) {std::cout <<this->g[i*nxyz+0]<<" ";} std::cout<<"\n";
			//std::cout<<"g 1 ";for (int i=0;i<ph_gcomp;i++) {std::cout <<this->g[i*nxyz+1]<<" ";} std::cout<<"\n";
			for (size_t i=0;i<gcompn.size();i++) {
				if (gcompn[i]=="H2O(g)") {this->iGwater = i;}
			}
		}
	std::cout << "c ";for (int i=0; i<ph_ncomp;i++){std::cout<< this->c[i*nxyz+0]<<" ";} std::cout<<"\n";
	std::cout <<"end phq init "<<"\n";

		/* Set boundary condition
		std::vector<double> bc_conc;
		std::vector<double> bc_f1; //bc_conc, 
		std::vector<int> bc1, bc2;
		int nbound = 2;
		bc1.resize(nbound, 0);                      // solution 0 from Initial IPhreeqc instance to boundary condition
		bc2.resize(nbound, 1);                     // no bc2 solution for mixing
		bc_f1.resize(nbound, 1.0);                  // mixing fraction for bc1
		status = phreeqc_rm->InitialPhreeqc2Concentrations(this->solu_conc, bc1, bc2, bc_f1); */
		
		//printing to file
		status = phreeqc_rm->CloseFiles();
		//status = phreeqc_rm->SetPrintChemistryOn(true, false, true); // workers, initial_phreeqc, utility

	}
	
	//--------------------------- RUN ------------------------------
	
	void run()
	{
		IRM_RESULT status;
		PhreeqcRM* phreeqc_rm = this->PhreeqcRM_ptr;
		
		status = phreeqc_rm->SetErrorHandlerMode(1);
		//status = phreeqc_rm->SetFilePrefix("phq_run");
		//phreeqc_rm->OpenFiles();
		//std::vector<int> print_chemistry_mask;
		//print_chemistry_mask.resize(nxyz, 1);//print_chemistry_mask[0]=1;print_chemistry_mask[50]=1;
		//status = phreeqc_rm->SetPrintChemistryMask(print_chemistry_mask); 
		//status = phreeqc_rm->SetPrintChemistryOn(true, false, false); // workers, initial_phreeqc, utility
		const std::vector<std::string> &components = phreeqc_rm->GetComponents();
		const std::vector<std::string> &gcompn = phreeqc_rm->GetGasComponents(); std::cout<<"phq ngcomp sz "<<gcompn.size()<<"ph ws "<<this->wsat[0]<<"\n";
		nxyz = this->data[0];
		std::cout << "in phq nxyz "<<nxyz<<" size c "<<this->c.size()<<" poro "<<this->poro[2]<<"\n";
		std::cout << "c 0 ";for (int i=0; i<ph_ncomp;i++){std::cout<< this->c[i*nxyz+2]<<" ";} std::cout<<"\n";
		status = phreeqc_rm->SetTemperature(this->temp);std::cout<<"temp done ";
		if (gcompn.size()>0) {
			std::cout<<"nb gcomp before "<<gcompn.size();//<<" gvol 3 "<<this->gvol[0]<<" sat "<<this->wsat[0]<<"\n";
			//std::cout<<"g 0 ";for (int i=0;i<ph_gcomp;i++) {std::cout <<this->g[i*nxyz+0]<<" ";} std::cout<<"\n";
			//std::cout<<"gm 0 ";for (int i=0;i<ph_gcomp;i++) {std::cout <<this->gm[i*nxyz+0]<<" ";} std::cout<<"\n";
			//status = phreeqc_rm->SetSaturation(this->wsat); // keep this order to set values
			status = phreeqc_rm->SetGasPhaseVolume(this->gvol);
			status = phreeqc_rm->SetGasCompMoles(this->gm);
			//status = phreeqc_rm->SetPressure(this->p); // does not seem to be usefull but we left it
			}
		else 
			{std::cout<<"phq wsat "<<this->wsat[0];
			//status = phreeqc_rm->SetSaturation(this->wsat);
			}
		std::cout<<"wsat done ";
		if (this->EK) {status = phreeqc_rm->SpeciesConcentrations2Module(this->spc);}
		else {status = phreeqc_rm->SetConcentrations(this->c);}
		std::cout<<"c done ";
		status = phreeqc_rm->SetTimeStep(this->tstep);std::cout<<"tstep "<<this->tstep<<"\n";
		//----- run here ------
		status = phreeqc_rm->RunCells();
		
		//status = phreeqc_rm->GetSaturation(this->wsat);
		//this->temp = phreeqc_rm->GetTemperature();
		if (this->EK) {status = phreeqc_rm->GetSpeciesConcentrations(this->c);}
		else {status = phreeqc_rm->GetConcentrations(this->c);}
	std::cout << "c 0 ";for (int i=0; i<ph_ncomp;i++){std::cout<< this->c[i*nxyz]<<" ";} std::cout<<"\n";
		if (gcompn.size()>0) {
			status = phreeqc_rm->GetGasCompPressures(this->g);
			status = phreeqc_rm->GetGasCompMoles(this->gm);
			status = phreeqc_rm->GetGasPhaseVolume(this->gvol);
			this->p = phreeqc_rm->GetPressure();
			}
		//writes to the dump file
		//bool dump_on = true;
		//bool append = false;
		//status = phreeqc_rm->SetDumpFileName("phreeqc.dmp");
		//status = phreeqc_rm->DumpModule(dump_on, append);
	//phreeqc_rm->SetSelectedOutputOn(true);
	std::cout << "c 0 ";for (int i=0; i<ph_ncomp;i++){std::cout<< this->c[i*nxyz+0]<<" ";} std::cout<<"\n";
	if (ph_gcomp>0) {
	//std::cout<<"nb gcomp after "<<gcompn.size()<<" gvol 0 "<<this->gvol[0]<<" sat "<<this->wsat[0]<<"\n";
		std::cout<<"g 0 ";for (int i=0;i<ph_gcomp;i++) {std::cout <<this->g[i*nxyz+0]<<" ";} std::cout<<"\n";
		std::cout<<"gm 0 ";for (int i=0;i<ph_gcomp;i++) {std::cout <<this->gm[i*nxyz+0]<<" ";} std::cout<<"\n";
		std::cout<<"p ";for (int i=0;i<5;i++) {std::cout<<this->p[i]<<" ";} std::cout<<"\n";
		}
	}
	
	void getSelOutput()
	{
		IRM_RESULT status;
		PhreeqcRM* phreeqc_rm = this->PhreeqcRM_ptr;
		int n_user = phreeqc_rm->GetNthSelectedOutputUserNumber(0);
		status = phreeqc_rm->SetCurrentSelectedOutputUserNumber(n_user);
		int ncol = phreeqc_rm->GetSelectedOutputColumnCount();std::cout<<"in selout ncol "<<ncol<<"\n";
		std::string heading;int nphase=0;
		for (int i = 0; i < ncol; i++)
			{
				status = phreeqc_rm->GetSelectedOutputHeading(i,heading);
				if (heading.substr(0,2)=="d_")  {nphase +=1;}
			}
		nxyz = this->data[0];std::cout<<"in selout nphase "<<nphase<<"\n";
		this->nselect = ncol-6-nphase;
		this->spc.resize((this->nselect)*nxyz,0);
		std::vector<double> so(ncol*nxyz);
		status = phreeqc_rm->GetSelectedOutput(so);
		//std::vector<std::string> headings;
		//headings.resize(ncol);
		for (int j = 0; j < nxyz; j++)
		{
			nphase = 0;
			for (int i = 6; i < ncol; i++)
			{
				status = phreeqc_rm->GetSelectedOutputHeading(i,heading);
				if (heading.substr(0,2)=="d_") {nphase +=1;}
				else 
					{
					this->spc[(i-6-nphase)*nxyz + j] = so[i*nxyz + j];
					if (j==10) {std::cout<< "selout sp "<< i << " : " << so[i*nxyz + j] << "\n";}
					}
			}
		}
	}
};
