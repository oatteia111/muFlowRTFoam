/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    muFlowRT

Description
    Stationary solver for incompressible single-phase flow in porous medium

Developers
    Olivier Atteia (the impes formulation come from Pierre Horgue in multiphaseporousFoam)

\*---------------------------------------------------------------------------*/
//#include <cmath.h>
//#include <iomanip> //NB when in < > don't add the .h
#include <stdlib.h>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <vector>
//#include <math.h>       /* atan */ don(t include it, it creates pb in opf when compiled with blue-cfd???
//#define MY_PI 3.14159265 // pb in blueCfd conflict with math?
#include <fstream>
#include <iterator>
#include <sstream>
#include <string>
#include <chrono>  // for high_resolution_clock
// #include "fvCFD.H" fro version 8-10, has disappeared now
#include "fvMesh.H"
#include "fvc.H"
#include "fvm.H"
#include "fvMatrices.H"
#include "IOdictionary.H"
#include "Time.H"
#include "argList.H"
#include "fvPatchFields.H"

#include "fvSchemes.H"
//#include "incompressiblePhase.H"
#include "fixedValueFvPatchFields.H"
#include "zeroGradientFvPatchFields.H"
#include "uniformDimensionedFields.H"
#include "inletOutletFvPatchField.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvPatchFieldMapper.H"
// #include "fvOptions.H"
#include "simpleControl.H"
#include "cellSet.H"

#include "phreeqc/initPhreeqc.H"
#include "partage.h"

// static allows to have ade finiotn ofr muFlow and another for phreeqc
#include <unistd.h> 
std::string get_current_dir() {
   char buff[FILENAME_MAX]; //create string buffer to hold path
   GetCurrentDir( buff, FILENAME_MAX );
   std::string current_working_dir(buff);
   return current_working_dir;
}
std::string cur_dir = get_current_dir();

std::vector<double> a(12,0.);

std::vector<double> c_ph,gm_ph,g_ph,poro,t_ph,foc_ph,p_ph,gvol,wsat,ractive,solu_conc,gas_conc,solu_species,Vmol,Ggrd;
std::vector<double> species;
std::vector<int> immobile;
std::vector<Foam::scalar> wTimes;
float atmPa=101325.;float pi=3.141592654;
float vmw,Cgtot,Gmtot,dure;
Foam::scalar dtForC,dtForChem,tnext,tunits;
int timeScheduleVersion = 1;
int i,j,iw,oindex,bindex,nsel;int rSteps=1;		     		  
int nxyz, ph_ncomp, ph_nspc, ph_gcomp, ph_nsolu; 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;//utilisteias and plugins declaration
#include "utilities.h" // for reading binary reading tables..
#include "myFunc.H"
// #include "transport/adaptiveReactiveDdtScheme.H" // new ddt solver to switch, bof added to matrix direclty

int main(int argc, char *argv[])
{
	my_phq freak;
	
	//init openFoam
	#include "foamVersion.H"
//#if defined(WM_PROJECT_VERSION_NUMBER) &&  WM_PROJECT_VERSION_NUMBER >= 10
  #define TIME_NAME(t) (t).name()
//#else
//  #define TIME_NAME(t) (t).timeName()
//#endif
	#include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    timeScheduleVersion = runTime.controlDict().lookupOrDefault<label>("timeScheduleVersion", 1);
    if (timeScheduleVersion != 1 && timeScheduleVersion != 2)
    {
        FatalErrorInFunction << "Unsupported timeScheduleVersion " << timeScheduleVersion
            << exit(FatalError);
    }
    if (!std::isfinite(runTime.value()) || runTime.value() < 0
        || !std::isfinite(runTime.endTime().value()) || runTime.endTime().value() <= 0
        || runTime.endTime().value() >= GREAT || runTime.value() > runTime.endTime().value())
    {
        FatalErrorInFunction << "Require 0 <= startTime <= endTime < GREAT, with positive endTime"
            << exit(FatalError);
    }
	//#include "polyMesh.H" was supposed to be usefull for cellLevles but creates pb elsewhere
    #include "readGravitationalAcceleration.H"
    #include "createFields.H"
    if (activateReaction == 1 && (!std::isfinite(reactionSteps)
        || reactionSteps == 0 || std::floor(reactionSteps) != reactionSteps
        || std::abs(reactionSteps) >= std::numeric_limits<int>::max()))
    {
        FatalErrorInFunction << "reactionSteps must be a nonzero integer; adaptive chemistry is not supported"
            << exit(FatalError);
    }
    if (activateReaction == 1 && runTime.value() != 0)
    {
        FatalErrorInFunction << "Reactive restart requires PHREEQC state restoration and is not supported"
            << exit(FatalError);
    }
	#include "readPicardControls.H"
	#include "flow/createThetaFields.H"
	//parms of the De=f(T) function
	//Info << "fDe parms "<< fDe_T.fparms[0] << " "<<fDe_T.fparms[1] << endl;
	#include "transport/createCFields.H"
	#include "transport/createTFields.H"
	#include "EK/createEKFields.H"
	
	#include "readCoupling.H"
	scalar deltaTFact = 1;
	// reading files times
	std::ifstream inputwTimes{cur_dir+"/constant/options/writetimes" };
	wTimes = {std::istream_iterator<scalar>{inputwTimes}, {}};
	if (wTimes.size() < 2 || !inputwTimes.eof())
	{
		FatalErrorInFunction
			<< "writetimes must contain a time-unit index and at least one output time"
			<< exit(FatalError);
	}
    if (wTimes[0] < 0 || wTimes[0] > 4 || std::floor(wTimes[0]) != wTimes[0])
    {
        FatalErrorInFunction << "Invalid time-unit index in writetimes" << exit(FatalError);
    }
    scalar lastWriteTime = 0;
    for (std::size_t index=1; index<wTimes.size(); ++index)
    {
        if (!std::isfinite(wTimes[index]) || wTimes[index] <= lastWriteTime || wTimes[index] >= GREAT)
        {
            FatalErrorInFunction << "Output times must be finite, positive, strictly increasing and below GREAT"
                << exit(FatalError);
        }
        lastWriteTime = wTimes[index];
    }
    if (timeScheduleVersion == 2 && wTimes.back() != runTime.endTime().value())
    {
        FatalErrorInFunction << "The final output must equal endTime" << exit(FatalError);
    }
	Info<<" wt0 "<<wTimes[0]<<endl;
		
	scalar tunits=1; float lunits=1.;
	std::vector<scalar> tu={1,60,3600,86400,3153600};tunits = tu[wTimes[0]]; // we need time units to send time to phrq in seconds and for pressure units
	wTimes[0] = 0; // required to calculate the reactive time sptes for the 1st time step 
	std::vector<float> lu={0.01,1,1000,0.325};lunits = lu[lg_units]; // we need length units for atmPa (lunit sis in transportproperties
	atmPa = atmPa*tunits*tunits/lunits;
	//cps = cps*tunits*tunits;cpw = cpw*tunits*tunits;lbdaTw = lbdaTw*tunits*tunits*tunits;lbdaTs = lbdaTs*tunits*tunits*tunits; // already mutilplied in input
	std::cout<<" tunits "<<tunits<<" atmpa "<<atmPa<<"\n";
	#include "flow/create2phaseFields.H"
	
#include "plugins/plugin_H.H" //variables to be modified before the H equation
#include "plugins/plugin_PS.H"
#include "plugins/plugin_Cgi.H"
	
	plugin_H plugH;
	plugin_PS plugPS;
	plugin_Cgi plugCgi;

	//########################## observations 
	//remove all files in the folder observation
	if (mesh.time().value()==0) {fname = cur_dir+"/observation"; if (fexists(fname)) {deleteFilesInDirectory(fname);} }
	//Observations : obspts a file with the name of each point and its x,y,z, coordinates
	fname = cur_dir+"/constant/options/obspts";int nobs = 0;std::vector<int> icello;outTable observ;
	if (fexists(fname))
		{
		observ = readTable(fname);nobs = observ.nrow;int ncol = observ.ncol;icello.resize(nobs);std::cout<<"nobs "<<nobs<<" "<<ncol<<"\n";
		for (int io=0;io<nobs;io++) {
			int ix=observ.data[io*ncol];int iz=observ.data[io*ncol+1];
			icello[io]=iz*ncell_lay+ix;std::cout<<"obs i "<<ix<<" "<<iz<<" "<<icello[io]<<"\n";
			}
		}	
	// observations : the files containing the number of the variable to be printed obsFlow, obsTrans, obsChem
	std::vector<scalar> obsFlow,obsTrans,obsSolu,obsGas;
	std::vector<int> obsFlowIndx,obsTransIndx,obsChemIndx;
	std::ifstream inputOflow{cur_dir+"/constant/options/obsFlow" }; 
	obsFlowIndx = {std::istream_iterator<int>{inputOflow}, {}};
	std::ifstream inputOtrans{cur_dir+"/constant/options/obsTrans" }; 
	obsTransIndx = {std::istream_iterator<int>{inputOtrans}, {}};
	std::ifstream inputOchem{cur_dir+"/constant/options/obsChem" }; //two values 0 or 1 for solutes and gases print
	obsChemIndx = {std::istream_iterator<int>{inputOchem}, {}}; 
	Info<<"end obs "<<endl;
	// budget variables
	std::vector<scalar> budFlow,budC,budT,budSolu,budGas;
	
	std::ifstream inputRactive{cur_dir+"/constant/options/ractive" }; // version 0 shall contain 0 for inactive and 1 for active reaction cell
	ractive = {std::istream_iterator<int>{inputRactive}, {}};
	std::vector<int>rinactive;
	std::ifstream inputInactiveCells{cur_dir+"/constant/options/inactiveCells" };
	inactiveCells = {std::istream_iterator<int>{inputInactiveCells}, {}};
	fname=cur_dir/"phqfoam.txt";std::ifstream inputData{fname};
	std::vector<int> ph_data{std::istream_iterator<int>{inputData}, {}}; //for (int i=0; i<7;i++){Info << "init nb "<< ph_data[i] << endl;}
	//--------------- for mpi -----------------
	//#include "add_mpi.h"

	//##############  phreeqc initialisation using initchem and data in phqfoam
	if (activateReaction==1)
{	
	rSteps = reactionSteps*1;
	//##################### make the initialization of phreeqc : data,poro, gvol
	nxyz=ph_data[0];ph_ncomp=ph_data[1];ph_gcomp=ph_data[2];ph_nsolu=ph_data[3]; //!!! nxyz here is inside the ractive part
	std::vector<int>temp(ncell,1);
	for (j=1;j<nxyz;j++) {temp[ractive[j]]=0;}
	for (j=1;j<ncell;j++) {if (temp[j]==1) {rinactive.push_back(j);} } 
	//Info<<"n cell "<<ncell<<" nxyz "<<nxyz<<" ract.size, ract(0) "<<ractive.size()<<" "<<ractive[0]<<" rinact.size, ract(0) "<<rinactive.size()<<" "<<rinactive[0]<<endl;
	if (activateEK) {freak.EK=true;} else {freak.EK=false;}
	std::cout<<cur_dir/"phreeqc.dat"<<"\n";
	std::string dbfile = cur_dir/"phreeqc.dat";
	freak.setDB(dbfile); //cur_dir/"phreeqc.dat");
	freak.setData(ph_data); //here we include the phqfoam data in freak it will be used by initphreeqc
	std::string chfile = cur_dir/"initChem.pqi";
	freak.setChemFile(chfile); //Info << "initCh read " << endl;
	//initiate poro and gas volume
	poro.resize(nxyz,0.);t_ph.resize(nxyz,0.);foc_ph.resize(nxyz,0.);p_ph.resize(nxyz,0.);Vmol.resize(nxyz,24);
	for (i=0;i<nxyz;i++) {poro[i]=eps[i];t_ph[i]=T[i];foc_ph[i]=foc[i];}
	wsat.resize(nxyz,1-1e-4);gvol.resize(nxyz,1);
	freak.setPoro(poro);freak.setWsat(wsat);freak.setTemp(t_ph);freak.setFoc(foc_ph);freak.setGvol(gvol);
	if (ph_gcomp>0) {
		// gases come from gases file
		std::ifstream inGases(cur_dir/"constant/options/gases");
		gas_conc = {std::istream_iterator<double>{inGases}, {}}; 
		inGases.close();

		gm_ph.resize(nxyz*ph_gcomp,0);
		//find the start index of gases in phqfoam file (ph_data variable)
		int start = 5;
		for (i=0;i<4;i++) {
			if (ph_data[start] == -1) {start += nxyz+1;} else {start += 2;}
			//Info<<"i data "<<i<<" start "<<start<<" value "<<ph_data[start]<<endl;
			}
		int igph = 0;
		Info<<"before phqinit"<<endl;
		for (int j=0;j<ncell;j++) {
			//p_ph[j]=p[j]/atmPa;Info<<"p_ph "<<p_ph[j];
			Vmol[j] = 24.46*(273.15+25.)/(273.15+T[j])/(p[j]/atmPa); 
			if (ph_data[start]==-1) {igph=ph_data[start+1+j];}
			for (int i=0;i<ph_gcomp;i++) {gm_ph[i*nxyz+j]=gas_conc[igph*ph_gcomp+i]/Vmol[j];Info<<" "<<igph<<" "<<gm_ph[i*nxyz+j];}
			Info<<endl;
			}
		freak.setGm(gm_ph);//freak.setP(p_ph);
		}
	//***first init of phreeqc
	int a0= freak.phqInit(); //if gas is present here the equil is not correct, it is fixed pressure(gas phase from phqfoam)
    if (a0 < 0)
    {
        FatalErrorInFunction << "PHREEQC initialization failed (status " << a0 << ")"
            << exit(FatalError);
    }
	a0=freak.getSelOutput();
    if (a0 < 0)
    {
        FatalErrorInFunction << "Cannot read PHREEQC initialization output" << exit(FatalError);
    }
	nsel = freak.nselect;//std::cout<<"1st phq, nsel "<<nsel<<" nxyz "<<nxyz<<"\n";
	species.resize(nxyz*nsel);
	for (size_t k;k<freak.spc.size();k++) {species[k]=freak.spc[k];} // put the starting concentrations

	
	//a0 = phqRun(freak);
	
	//################ writes the solutions to file + compnames
	//solutions are obtained from boundary conditions
	//these solutions will be used for chem BCs
	std::ofstream outSolu(cur_dir/"constant/options/solutions");
	solu_conc.resize(ph_nsolu*ph_ncomp,0.);Info << "nsolu "<<ph_nsolu << " ncomp "<< ph_ncomp <<" bcs "<<freak.bc_conc.size()<<endl;
	//if (ph_gcomp==0) { // seems to work only for solutions
	for (i=0;i<ph_nsolu;i++) // solu number
		{ 
		for (j=0;j<ph_ncomp;j++) // component number
			{
			float a = freak.bc_conc[j*ph_nsolu+i];std::cout<<a<<"\n";
			solu_conc[i*ph_ncomp+j] = a; outSolu << a << "\n"; 
			} 
		}		
	outSolu.close();
	// if activate EK stores the species
	if (activateEK)
		{
		solu_species.resize(ph_nsolu*freak.nspc,0.);Info << "nsolu "<<ph_nsolu << " nspecies "<< freak.nspc <<endl;
		for (i=0;i<ph_nsolu;i++) // solu number
			{ 
			for (j=0;j<freak.nspc;j++) // component number
				{
				float a = freak.bc_species[j*ph_nsolu+i];
				solu_species[i*freak.nspc+j] = a; //outSolu << a << "\n"; 
				} 
			}		
		}
	//##############" build the c_ph and gm_ph fields and get conc from phreeqc (c_ph=freak.c but needed two variables for format questions)
	// get conc calculated by phreeqc from initchem and phqfoam
	c_ph.resize(nxyz*ph_ncomp,0);
	for (i=0; i<ph_ncomp;i++)
		for (int j=0;j<nxyz;j++) {c_ph[i*nxyz+j] = freak.c[i*nxyz+j];}

	gm_ph.resize(nxyz*ph_gcomp,1e-9); // moles of gas
	//in order to have correct correct pressure, we set gm by considering gvol
	std::cout<<" ng comp "<<ph_gcomp<<"\n";
	if (ph_gcomp>0)
		{
		std::cout<<" press start "<<p[0]/atmPa<<"\n";;
		for (j=0;j<nxyz;j++) 
			{
			gvol[j] = max(eps[j]*(1-sw[j]),1e-4); //
			Gmtot = 0;
			//for (i=0;i<ph_gcomp;i++) {Gmtot += freak.gm[i*nxyz+j];}
			if (j<15) {std::cout<<"p "<<p[j]/atmPa<<"gv "<<gvol[j];} //<<" Vm "<<Vmol[j];
			for (i=0;i<ph_gcomp;i++) {gm_ph[i*nxyz+j] = freak.gm[i*nxyz+j]*gvol[j];  //*p[j]/atmPa;
				if (j<15) {std::cout<<" "<<gm_ph[i*nxyz+j];}
				} //phreeqc sends back moles for 1L RV
			if (j<15) {Info<<endl;}
			//freak.gm[i*nxyz+j]/Gmtot*gvol[j]/Vmol[j]			
			}
		//just print c
		Info<<"conc"<<endl;
		for (j=0;j<nxyz;j++) 
			{
			for (i=0;i<ph_ncomp;i++) {if (j<15) {std::cout<<" "<<freak.c[i*nxyz+j];}} 
			if (j<15) {Info<<endl;}
			}
			
		std::cout<<"end get gm \n";
		freak.setGvol(gvol); // set gas volume in phreeqc
		freak.setGm(gm_ph);//transfer gm_ph to freak
		freak.setTemp(t_ph);
		p_ph.resize(nxyz);
		for (j=0;j<nxyz;j++) {p_ph[j]=p[j]/atmPa;}
		freak.setP(p_ph); //not possible to set pressure and volume
		a0= freak.phqRun(); //****PHQ RUN with equilibration with true gas phase
        if (a0 < 0)
        {
            FatalErrorInFunction << "Initial gas equilibration failed (status " << a0 << ")"
                << exit(FatalError);
        }
		//(recalculate Vm) no, just to print
		for (j=0;j<nxyz;j++)  {
			//Gmtot = 0;
			//for (i=0;i<ph_gcomp;i++) {Gmtot += freak.gm[i*nxyz+j];}
			//Vmol[j] = gvol[j]/Gmtot; //
			std::cout<<"p "<<freak.p[j]<<"gv "<<freak.gvol[j];
			for (i=0;i<ph_gcomp;i++) {Info<<" "<<freak.gm[i*nxyz+j];} Info<<endl;
			}
		}
	
	iw = freak.iGwater;	
	// set values of cells when restart
	
	
	Info<<"end c_ph and gm_ph "<<endl;
	
	//###################"" loading immobile component if present  ###############""
	std::vector<std::string> immobStr;
	std::ifstream inputImmobile{cur_dir+"/constant/options/immobile" }; // version 0 shall contain 0 for inactive and 1 for active reaction cell
	immobStr = {std::istream_iterator<std::string>{inputImmobile}, {}};
	immobile.resize(ph_ncomp,0);
	for (i=0;i<immobStr.size();i++) 
		{
		int ic = std::find(freak.comp.begin(),freak.comp.end(),immobStr[i])-freak.comp.begin(); //find the position of str in comp
		immobile[ic]=1;std::cout<<"immob "<<immobStr[i]<<" "<<ic<<" "<<immobile[ic]<<"\n";
		}
	Info<<"end immobile "<<endl;
	
} //end of activateReation
	
	else  //only flow or flow+transport
	{
		ph_ncomp=0;ph_gcomp=0;freak.gcomp={""}; // this is needed to initate pluginCgi (not very clean)
	}
	Info <<"ncomp "<<ph_ncomp<<endl;

	#include "transport/createCwiFields.H"
	#include "transport/createCgiFields.H"
	// #include "flow/update2phaseFields.H"									 
	
	//######################## run the steady state for hp
	dimensionedScalar st = runTime.startTime();
	//dimensionedScalar et = runTime.endTime();
	scalar dt0 = mesh.time().deltaTValue();
	scalar residu;//float dt0=wtime/50.;Info<<"wtime 0 "<<wtime<<" dt "<<dt0<<endl;
	if ((flowStartSteady==1)&&(flowType>0)&&(flowType<=2))
		{
		runTime.setDeltaT(dt0);
		#include "flow/hstdEqn.H"
		}
	// ###################  starting timer ######################
	Info <<"st time "<<st<<endl;
	//runTime.setEndTime(et); Info<<"end "<<runTime.endTime()<<endl;
	runTime.runTimeModifiable();
	Info <<"dt "<<dt0<<endl;
	
	runTime.setDeltaTNoAdjust(dt0);
	time=mesh.time().value();
	int tstep=runTime.timeIndex(),tcnt=0,rcnt=1,iflowStep=0;
	std::size_t itwstep = 1; // the first entry stores the time unit
	int flagW=0;
	scalar presentTime=runTime.value(),oldTimeReac=runTime.value();
	scalar oldTimeTransp=runTime.value(),nextTimeTransp=runTime.value(),newDeltaT=dt0;
	// Restart from the first output after the state already loaded at startTime.
	while (itwstep < wTimes.size() && wTimes[itwstep] <= runTime.value()) {++itwstep;}
	scalar wtime = itwstep < wTimes.size() ? wTimes[itwstep] : GREAT;
	//restart for solid species in phreeqc (get minerals from species file, equilibrate surface and exchange)
	/*
	if ((itwstep>1)&&(activateReaction==1)) 
		{
		std::vector<double> g_ph(ph_gcomp*ncell);										   
		//# include "transport/createCwiFields.H"
		forAll(Cw,i) {for (j=0; j<nxyz;j++) {c_ph[i*nxyz+j] = Cw[i]()[ractive[j]];} };std::cout<<"conc read "<<Cw[4]()[233]<<"\n";
		if (ph_gcomp>0)
			 forAll(Cg,i) {for (j=0; j<nxyz;j++) {g_ph[i*nxyz+j] = Cg[i]()[ractive[j]];} };//std::cout<<"conc read "<<Cw[4]()[233]<<"\n";				 																														
		int a0 = freak.phqRestart(ph_data,std::to_string(int(time)),c_ph,g_ph); // I did not find a way to send Cw -> dimensoin error?
		fname=cur_dir/"phqfoam1.txt";std::ifstream inputData1{fname};
		std::vector<int> ph_data{std::istream_iterator<int>{inputData1}, {}}; //for (int i=0; i<7;i++){Info << "init nb "<< ph_data[i] << endl;}
		freak.setData(ph_data);
		freak.setChemFile(cur_dir/"initChem1.pqi"); //Info << "initCh read " << endl;
		a0=freak.phqInit();
		forAll(Cw,i) {for (j=0; j<nxyz;j++) {Cw[i]()[ractive[j]]=freak.c[i*nxyz+j];} };std::cout<<"conc read "<<Cw[4]()[233]<<"\n";
		itwstep+=1;
		std::cout<<"end restart "<< Cw[4]()[0]<<" \n";
		}
		*/
	Info<<"wtime "<<itwstep<<" "<<wtime<<endl;
	Info<<"time rebuilt st "<<runTime.startTime()<<" dt "<<runTime.deltaTValue()<<endl;
	
	//################### attribute to Cwi and Cgi concentrations/pressures from phreeqc
	int icnt = 0;
	//Info<<"ractive 2 "<<ractive[2]<<endl;
	//Info<<" gmph 2 "<<  gm_ph[2]<<" frk.c size "<<freak.c.size()<<endl;
	
	if (activateReaction==1) 
	{
		forAll(Cw,i) {
			for (j=0; j<nxyz;j++){Cw[i]()[ractive[j]] = freak.c[i*nxyz+j];} // transfer freak.c to Cw
			for (j=1;j<ncell-nxyz;j++) {Cw[i]()[rinactive[j]]=solu_conc[i];}
			} 
		if (ph_gcomp>0) { // there are gases only when reaction is present, the freak.g transmit pressures in bars
			Info<<"assign gases to Cg"<<endl;
			for (j=0; j<nxyz;j++)
				{
				//Cgtot = 0;
				//forAll(Cg,i){Cgtot += freak.g[i*nxyz+j];} 
				for (i=0;i<ph_gcomp;i++) {Cg[i]()[j] = freak.gm[i*nxyz+j]/gvol[j];} // phreeqc is in mol/RV
				if (j<15) {for (i=0;i<ph_gcomp;i++) {Info<<Cg[i]()[j]<<" ";} Info<<endl;}
				}
		}
	}
	//if (ph_gcomp>1) {Info<<" cg 0 1 "<<Cg[0]()[1]<<endl;}
	//Info<<"Cg "; for (j=0;j<nxyz;j++)  {for (i=0;i<ph_gcomp;i++) {Info<<" "<<Cg[i]()[j];} Info<<endl;}

	
	//--------------------------------------------------------------------------------
	//                          MAIN LOOP
	//---------------------------------------------------------------------------------
	
	//###########################  plugins
	plugH.init(cur_dir,transportProperties,mesh,runTime,listCouples,lunits,tunits);
	plugPS.init(cur_dir,transportProperties,mesh,freak,listCouples);
	plugCgi.init(cur_dir,transportProperties,mesh,freak,listCouples);
    while (true)
    {
        presentTime = runTime.value();
        if (presentTime < runTime.endTime().value())
        {
            #include "advanceBoundaries.H"
        }
        const scalar writeTime = wtime;
        const bool hasWriteTime = itwstep < wTimes.size();
        const std::size_t intervalIndex = std::min(itwstep, wTimes.size()-1);
        const scalar intervalStart = wTimes[intervalIndex-1];
        const scalar reactStep = activateReaction == 1 && rSteps > 0
            ? (wTimes[intervalIndex]-intervalStart)/rSteps : scalar(0);
        const scalar reactionTime = reactStep > 0 && hasWriteTime
            ? (rcnt == rSteps ? writeTime : intervalStart+reactStep*rcnt) : scalar(GREAT);
        scalar stepTime = presentTime;
        const auto prepareStep = [&]()
        {
            const scalar endTime = runTime.endTime().value();
            if (presentTime >= endTime) {return;}
            maxDeltaT = runTime.controlDict().lookupOrDefault<scalar>("maxDeltaT", GREAT);
            minDeltaT = runTime.controlDict().lookupOrDefault<scalar>("minDeltaT", 1e-8);
            if (!std::isfinite(newDeltaT) || newDeltaT <= 0
                || !std::isfinite(minDeltaT) || !std::isfinite(maxDeltaT)
                || minDeltaT <= 0 || maxDeltaT < minDeltaT)
            {
                FatalErrorInFunction << "Invalid numerical timestep limits at " << presentTime
                    << exit(FatalError);
            }
            try
            {
                stepTime = muFlowTime::stepEnd
                (
                    presentTime, min(max(newDeltaT,minDeltaT),maxDeltaT),
                    writeTime, tnext, reactionTime, endTime
                );
            }
            catch (const std::exception& error)
            {
                FatalErrorInFunction << error.what() << " at " << presentTime << exit(FatalError);
            }
            runTime.setDeltaTNoAdjust(stepTime-presentTime);
        };
        // run() uses deltaT in its stop test; clip before it, then refresh controls.
        prepareStep();
        if (!runTime.run()) {break;}
        prepareStep();
        if (presentTime >= runTime.endTime().value()) {continue;}
        if (timeScheduleVersion == 2 && wTimes.back() != runTime.endTime().value())
        {
            FatalErrorInFunction << "Changing endTime requires regenerating writetimes"
                << exit(FatalError);
        }
        flagW = hasWriteTime && stepTime == writeTime;
        runTime++;tstep++;tcnt++;
        scalar dt = runTime.deltaTValue();
        Info << "time = " << runTime.value() << " deltaT = " << dt
            << " next boundary " << tnext << endl;

		//***********************  solve coupling case (we assume all coupling require h/C loop) *******************************
		if (coupling ==1) // we assume Picard
			{
			iterPicard = 0; resPicard = 1000.;float residu0 = 1000;deltaTFact = 1.;
			h1 = h;
			while (resPicard > tolPicard)
				{
				iterPicard++;
				h = h_tmp; 
				#include "flow/hEqn.H"
				if (activateReaction==0){
					#include "transport/CEqn.H"
				}
				else {
					forAll(Cw,i) {Cw[i]().storePrevIter();} // for cells outside calculation
					#include "transport/CwiEqn.H"
				}
				resPicard = gMax((mag(h-h_tmp))->internalField());
				if (residu>residu0) {h  = h*0.75 + h.prevIter()*0.25;}//relax nb residu is calc in hEqn
				h_tmp = h;
				residu0 = residu;
				if (iterPicard == maxIterPicard)
					{ break; }
				}  // end picard iter
			Info << "Picard nb iterations : "<<iterPicard<<endl;
			if (activateReaction==1 && resPicard <= tolPicard) {
				#include "phreeqc/calcReaction.H"
			}
			if (!std::isfinite(resPicard) || resPicard > tolPicard)
            {
                FatalErrorInFunction << "Coupled Picard solve failed after accepted time "
                    << presentTime << "; rejected step was not written" << exit(FatalError);
            }
			else if (iterPicard <= nIterStability) {deltaTFact = 1.2;}
			newDeltaT = min(deltaTFact*newDeltaT,maxDeltaT); //runTime.deltaTValue()
			if (iterPicard>1) {
				#include "flow/hEqn.H"
				}
			}

		//***********************  solve transient flow   *******************************
		if ((coupling ==0)&&(flowType>0) )
			{
			iflowStep++;
			#include "flow.H"
			if (ph_gcomp>0) 
				{
				gvol.resize(nxyz);
				for (j=0; j<nxyz;j++) {gvol[j]=max(eps[j]*(1-sw[j]),1e-4);} 
				}
			}
		
		//***************  solve Transport & reaction (no coupling) *************************
		//float sumT=0;float sumPF=0;float sumPFa=0;

		if ((coupling ==0)&&(activateThermal==1))
			{
			#include "transport/TEqn.H"
			}
					
		if ((coupling ==0)&&(activateEK==1)) 
			{
			#include "EK/NernstPlanck.H"
			}
		
		if ((coupling ==0)&&(activateEK==0) && (activateTransport==1)) {
			//if ((mesh.time().value()>=nextTimeTransp)||(iflowStep>20)) {
				if (activateReaction==0) {
					#include "transport/CEqn.H"
					}
				else { //reaction occurs
					forAll(Cw,i) {Cw[i]().storePrevIter();} // for cells outside calculation
					#include "transport/CwiEqn.H"
					if (ph_gcomp>0) {
						forAll(Cg,i) {Cg[i]().storePrevIter();}
						#include "transport/CgiEqn.H"
						}
					#include "phreeqc/calcReaction.H"
					}
				/* for variable transp steps
				iflowStep = 0;
				nextTimeTransp = mesh.time().value()+min(dtForC/2.,maxDeltaT);
				oldTimeTransp = runTime.value()*1;
				std::cout<<"end trsp, pres "<<mesh.time().value()<<" next "<<nextTimeTransp<<"\n";
				*/
				}
			//}
			
		// Use the requested time for output metadata, not accumulated step roundoff.
		runTime.setTime(stepTime,runTime.timeIndex());
		//bool ts = runTime.write();
		//write
		#include "observation.H"
		#include "budget.H"
		
		if (flagW==1) {
            if (!runTime.writeNow())
            {
                FatalErrorInFunction << "Failed to write output at " << writeTime << exit(FatalError);
            }
            Info<<"writing "<<writeTime<<endl;
			if (rSteps>0) {tcnt=0;rcnt=1;}
			}
		
		//if (flowType==4) {phiGr.write();}
		if (activateReaction==1  && flagW==1) {
			phiw.write();phig.write();
			std::ofstream outFile(cur_dir/TIME_NAME(mesh.time())/"Species");
			outFile.unsetf(std::ios::scientific);outFile.precision(6);
			std::cout<<"write nsel "<<nsel<<" nxyz "<<nxyz<<"\n";
			for (j=0;j<nxyz;j++)
				{ for (i=0;i<nsel;i++) {outFile << species[i*nxyz + j]<<" ";} outFile <<"\n"; }
            outFile.close();
            if (!outFile)
            {
                FatalErrorInFunction << "Failed to write Species at " << writeTime << exit(FatalError);
            }
			}

		if (flagW==1)
		{
			// Consume the target only after a successful step and its output.
			++itwstep;
			wtime = itwstep < wTimes.size() ? wTimes[itwstep] : GREAT;
			flagW=0;
		}
		Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
			<< "  ClockTime = " << runTime.elapsedClockTime() << " s"
			<< nl << endl;

		Info<< "tnext "<<tnext<<" End tstep\n" << endl;

	}
	
	Info<<"Normal termination of OpenFoam"<<endl;
	/*
	#ifdef USE_MPI
		freak.PhreeqcRM_ptr->MpiWorkerBreak();
		int status = MPI_Finalize();
	#endif
	*/
    return 0;
}

// ************************************************************************* //
