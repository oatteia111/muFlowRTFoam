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
#include <vector>
//#include <math.h>       /* atan */ don(t include it, it creates pb in opf when compiled with blue-cfd???
//#define MY_PI 3.14159265 // pb in blueCfd conflict with math?
#include <fstream>
#include <iterator>
#include <sstream>
#include <string>
#include <chrono>  // for high_resolution_clock
#include "fvCFD.H"
#include "fvSchemes.H"
//#include "incompressiblePhase.H"
#include "inletOutletFvPatchField.H"
#include "volFields.H"
#include "fvPatchFieldMapper.H"
// #include "fvOptions.H"
#include "simpleControl.H"
#include "cellSet.H"

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

std::vector<double> a(12,0.);
#include "phreeqc/initPhreeqc.H"
std::vector<double> c_ph,gm_ph,g_ph,poro,t_ph,foc_ph,p_ph,gvol,wsat,ractive,solu_conc,gas_conc,solu_species,Vmol,Ggrd;
std::vector<int> immobile;
std::vector<float> wTimes;
float atmPa=101325.;float pi=3.141592654;
float vmw,Cgtot,Gmtot,dtForC,dtForChem,tnext;
int i,j,iw,oindex,bindex;int rSteps=1;		     		  

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;//utilisteias and plugins declaration
#include "utilities.h" // for reading binary reading tables..
#include "myFunc.H"
#include "plugins/plugin_H.H" //variables to be modified before the H equation
#include "plugins/plugin_PS.H"
#include "plugins/plugin_Cgi.H"
// #include "transport/adaptiveReactiveDdtScheme.H" // new ddt solver to switch, bof added to matrix direclty

int main(int argc, char *argv[])
{
	my_phq freak; 
	plugin_H plugH;
	plugin_PS plugPS;
	plugin_Cgi plugCgi;

	//init openFoam
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
	//#include "polyMesh.H" was supposed to be usefull for cellLevles but creates pb elsewhere
    #include "readGravitationalAcceleration.H"
    #include "createFields.H"
	#include "readPicardControls.H"
	#include "flow/createThetaFields.H"
	//parms of the De=f(T) function
	//Info << "fDe parms "<< fDe_T.fparms[0] << " "<<fDe_T.fparms[1] << endl;
	#include "transport/createCFields.H"
	#include "transport/createTFields.H"
	#include "EK/createEKFields.H"
	
	#include "readCoupling.H"
	
	// reading files times
	std::ifstream inputwTimes{cur_dir+"/constant/options/writetimes" };
	wTimes = {std::istream_iterator<float>{inputwTimes}, {}};Info<<" wt0 "<<wTimes[0]<<endl;
	int tunits=1;float lunits=1.;
	std::vector<int> tu={1,60,3600,86400,3153600};tunits=tu[wTimes[0]]; // we need time units to send time to phrq in seconds and for pressure units
	wTimes[0] = 0; // required to calculate the reactive time sptes for the 1st time step 
	std::vector<float> lu={0.01,1,1000,0.325};lunits = lu[lg_units]; // we need length units for atmPa (lunit sis in transportproperties
	atmPa = atmPa*tunits*tunits/lunits;
	//cps = cps*tunits*tunits;cpw = cpw*tunits*tunits;lbdaTw = lbdaTw*tunits*tunits*tunits;lbdaTs = lbdaTs*tunits*tunits*tunits; // already mutilplied in input
	std::cout<<" tunits "<<tunits<<" atmpa "<<atmPa<<"\n";
	#include "flow/create2phaseFields.H"
	
	//########################## observations 
	//remove all files in the folder observation
	if (mesh.time().value()==0) {fname = cur_dir+"/observation"; if (fexists(fname)) {deleteFilesInDirectory(fname);} }
	//Observations : obspts a file with the name of each point and its x,y,z, coordinates
	fname = cur_dir+"/constant/options/obspts";int nobs;std::vector<int> icello;outTable observ;
	if (fexists(fname))
		{
		observ = readTable(fname);nobs = observ.nrow;int ncol = observ.ncol;icello.resize(nobs);std::cout<<"nobs "<<nobs<<" "<<ncol<<"\n";
		for (int io=0;io<nobs;io++) {
			int ix=observ.data[io*ncol];int iz=observ.data[io*ncol+1];
			icello[io]=iz*ncell_lay+ix;std::cout<<"obs i "<<ix<<" "<<iz<<" "<<icello[io]<<"\n";
			}
		}	
	// observations : the files containing the number of the variable to be printed obsFlow, obsTrans, obsChem
	std::vector<float> obsFlow,obsTrans,obsSolu,obsGas;
	std::vector<int> obsFlowIndx,obsTransIndx,obsChemIndx;
	std::ifstream inputOflow{cur_dir+"/constant/options/obsFlow" }; 
	obsFlowIndx = {std::istream_iterator<int>{inputOflow}, {}};
	std::ifstream inputOtrans{cur_dir+"/constant/options/obsTrans" }; 
	obsTransIndx = {std::istream_iterator<int>{inputOtrans}, {}};
	std::ifstream inputOchem{cur_dir+"/constant/options/obsChem" }; //two values 0 or 1 for solutes and gases print
	obsChemIndx = {std::istream_iterator<int>{inputOchem}, {}}; 
	Info<<"end obs "<<endl;
	// budget variables
	std::vector<float> budFlow,budC,budT,budSolu,budGas;
	
	std::ifstream inputRactive{cur_dir+"/constant/options/ractive" }; // version 0 shall contain 0 for inactive and 1 for active reaction cell
	ractive = {std::istream_iterator<int>{inputRactive}, {}};
	std::vector<int>rinactive;
	std::ifstream inputInactiveCells{cur_dir+"/constant/options/inactiveCells" }; // version 0 shall contain 0 for inactive and 1 for active reaction cell
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
	freak.setDB(cur_dir/"phreeqc.dat");
	freak.setData(ph_data); //here we include the phqfoam data in freak it will be used by initphreeqc
	freak.setChemFile(cur_dir/"initChem.pqi"); //Info << "initCh read " << endl;
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
			Info<<"i data "<<i<<" start "<<start<<" value "<<ph_data[start]<<endl;}
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
	int a0= phqInit(freak); //if gas is present here the equil is not correct, it is fixed pressure(gas phase from phqfoam)
	//a0 = phqRun(freak);
	
	//################ writes the solutions to file + compnames
	//solutions are obtained from boundary conditions
	//these solutions will be used for chem BCs
	std::ofstream outSolu(cur_dir/"constant/options/solutions");
	solu_conc.resize(ph_nsolu*ph_ncomp,0.);Info << "nsolu "<<ph_nsolu << " ncomp "<< ph_ncomp <<endl;
	//if (ph_gcomp==0) { // seems to work only for solutions
		for (i=0;i<ph_nsolu;i++) // solu number
			{ 
			for (j=0;j<ph_ncomp;j++) // component number
				{
				float a = freak.bc_conc[j*ph_nsolu+i];
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
		a0= phqRun(freak); //****PHQ RUN with equilibration with true gas phase
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
	dimensionedScalar et = runTime.endTime();
	float dt0 = mesh.time().deltaTValue();
	scalar residu;//float dt0=wtime/50.;Info<<"wtime 0 "<<wtime<<" dt "<<dt0<<endl;
	if ((flowStartSteady==1)&&(flowType>0)&&(flowType<=2))
		{
		runTime.setDeltaT(dt0);
		#include "flow/hstdEqn.H"
		}
	// ###################  starting timer ######################
	Info <<"st time "<<st<<endl;
	runTime.setEndTime(et); Info<<"end "<<runTime.endTime()<<endl;
	runTime.runTimeModifiable();
	Info <<"dt "<<dt0<<endl;
	
	runTime.setDeltaT(dt0);runTime.setTime(st,0);
	time=mesh.time().value();
	if (time<dt0) {runTime.setTime(dt0,0);} // not to have time=0, it hsould be  dt0
	int tstep=0;int itwstep = 1;int tcnt = 0;int rcnt = 0;float wtime=wTimes[itwstep];int flagW=0; // first tstep is 1 (first vlaue is tunits)
	float oldTime,oldTimeReac;float newDeltaT=dt0;float dteps=(wTimes[2]-wTimes[1])/1e5;int rewind=0;
	//search time step for restart (if start time is higher than first write time
	while (time>wTimes[itwstep]) {itwstep+=1;} 
	//restart for solid species in phreeqc (get minerals from species file, equilibrate surface and exchange)
	if ((itwstep>1)&&(activateReaction==1)) 
		{
		std::vector<double> g_ph(ph_gcomp*ncell);										   
		//# include "transport/createCwiFields.H"
		forAll(Cw,i) {for (j=0; j<nxyz;j++) {c_ph[i*nxyz+j] = Cw[i]()[ractive[j]];} };std::cout<<"conc read "<<Cw[4]()[233]<<"\n";
		if (ph_gcomp>0)
			 forAll(Cg,i) {for (j=0; j<nxyz;j++) {g_ph[i*nxyz+j] = Cg[i]()[ractive[j]];} };//std::cout<<"conc read "<<Cw[4]()[233]<<"\n";				 																														
		int a0 = phqRestart(freak, ph_data,std::to_string(int(time)),c_ph,g_ph); // I did not find a way to send Cw -> dimensoin error?
		fname=cur_dir/"phqfoam1.txt";std::ifstream inputData1{fname};
		std::vector<int> ph_data{std::istream_iterator<int>{inputData1}, {}}; //for (int i=0; i<7;i++){Info << "init nb "<< ph_data[i] << endl;}
		freak.setData(ph_data);
		freak.setChemFile(cur_dir/"initChem1.pqi"); //Info << "initCh read " << endl;
		a0=phqInit(freak);
		forAll(Cw,i) {for (j=0; j<nxyz;j++) {Cw[i]()[ractive[j]]=freak.c[i*nxyz+j];} };std::cout<<"conc read "<<Cw[4]()[233]<<"\n";
		itwstep+=1;
		std::cout<<"end restart "<< Cw[4]()[0]<<" \n";
		}
	wtime=wTimes[itwstep];Info<<"wtime "<<itwstep<<" "<<wTimes[itwstep]<<endl;
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
				forAll(Cg,i){Cg[i]()[j] = freak.gm[i*nxyz+j]/gvol[j];} // phreeqc is in mol/RV
				if (j<15) {forAll(Cg,i) {Info<<Cg[i]()[j]<<" ";} Info<<endl;}
				}
		}
	}
	//if (ph_gcomp>1) {Info<<" cg 0 1 "<<Cg[0]()[1]<<endl;}
	//Info<<"Cg "; for (j=0;j<nxyz;j++)  {for (i=0;i<ph_gcomp;i++) {Info<<" "<<Cg[i]()[j];} Info<<endl;}

	
	//--------------------------------------------------------------------------------
	//                          MAIN LOOP
	//---------------------------------------------------------------------------------
	
	//###########################  plugins
	plugH.init(cur_dir,transportProperties,mesh,runTime,listCouples);
	plugPS.init(cur_dir,transportProperties,mesh,freak,listCouples);
	plugCgi.init(cur_dir,transportProperties,mesh,freak,listCouples);
	int flagDeltaT;//newDeltaT = minDeltaT*10;
	scalar reactStep = (wTimes[1]-wTimes[0])/rSteps; // length of the reaction step
	
	while (runTime.run())
    {
		//****************** set time step to stop at writeTimes
		if (rewind==1) {runTime.setTime(oldTime,tstep);newDeltaT /=10.;rewind=0;} //rewind is when phreeqc makes error
		oldTime = mesh.time().value();
		int flg = 0;Info<<"dts : min "<<minDeltaT<<" max "<<maxDeltaT<<" new "<<newDeltaT<<endl;
		newDeltaT= min(max(newDeltaT,minDeltaT),maxDeltaT);
		/*
		if ((flagDeltaT==1)&&(itwstep<=wTimes.size())) {  // the previous time step showed that we reach now a change in BC condition
			newDeltaT = min(newDeltaT/20,(wTimes[itwstep+1]-wTimes[itwstep])/100); // for unsat wtimes can be quite distant
			flagDeltaT=0;
			}
		*/
		float dt1 = wtime - oldTime;  //to catch the writing time
		float dt2 = tnext - oldTime; //to catch the time when BC change
		Info<<" dt1 "<<dt1<<" dt2 "<<dt2<< "newdt "<<newDeltaT+dteps<<endl;
		if ((dt1 <= dt2)&&(dt1<=newDeltaT+dteps)&&(dt1>0)) //write
			{runTime.setDeltaTNoAdjust(dt1);itwstep+=1;wtime=wTimes[itwstep];
			flagW=1;flg=1;dteps=(wTimes[itwstep+1]-wTimes[itwstep])/1e4;
			newDeltaT /= 2.; tnext=runTime.endTime().value();
			}
		else if ((dt2 <= dt1)&&(dt2<= newDeltaT+dteps)&&(dt2>0)) //BC change
			{runTime.setDeltaTNoAdjust(dt2);tnext=runTime.endTime().value();
			newDeltaT = min(newDeltaT/20,(wTimes[itwstep+1]-wTimes[itwstep])/100);
			//flagDeltaT=1;
			flg=1;} //;tnext=runTime.endTime().value()
		//Info<<" flg "<<flg<<endl;
		else {runTime.setDeltaT(newDeltaT);} // (flg==0)  classical case
		if (dt1==0) {itwstep+=1;wtime=wTimes[itwstep];flagW=1;}
		//Info <<"newDeltaT "<<newDeltaT<<endl;
		Info<<"i time "<<itwstep<<" oldt "<<oldTime<<" wt "<<wtime<<" deltaT "<<float(runTime.deltaTValue())<<" flgW "<<flagW<<endl;
		runTime.read();
		runTime++;tstep++;
		float dt = runTime.deltaTValue();
		scalar reactStep = (wTimes[itwstep]-wTimes[itwstep-1])/rSteps; // length of the reaction step
		//if (mesh.time().value()==wtime) {flagW=1;}
		Info <<"time = "<< mesh.time().value() <<" deltaT = " <<  dt << " tnext "<<tnext<<" newdeltaT "<<newDeltaT<<" reactStep "<<reactStep<<endl;
		
		//***********************  solve transient flow   *******************************
		if (flowType>0) {
			#include "flow.H"
			}
		if (rewind==1) {continue;} //for case picard has failed on flow
		//Info<<"p after flow ";for (int j=0; j<nxyz;j++) {Info<<p[j]/atmPa<<" ";};Info<<endl;
		//Info<<"sw after flow ";for (int j=0; j<nxyz;j++) {Info<<sw[j]<<" ";};Info<<endl;
		if (ph_gcomp>0) {
			gvol.resize(nxyz);
			for (j=0; j<nxyz;j++) {gvol[j]=max(eps[j]*(1-sw[j]),1e-4);} 
			}
		
		//***************  solve Transport  *************************
		//float sumT=0;float sumPF=0;float sumPFa=0;

		if (activateThermal==1) {
			#include "transport/TEqn.H"
			}
					
		if (activateEK==1) {
			#include "EK/NernstPlanck.H"
			}
		else
		{
		if (activateTransport==1) {
			if (activateReaction==0) {
				#include "transport/CEqn.H"
				}
			
			else {
				forAll(Cw,i) {Cw[i]().storePrevIter();} // for cells outside calculation
				//in case of ractive reset all aout of reactive to iniital
				//Info<<"Cw mid "<<Cw[7]()[int(ncell/2)]<<endl;
				#include "transport/CwiEqn.H"
				//Info<<"Cw mid "<<Cw[7]()[int(ncell/2)]<<endl;

				if (ph_gcomp>0) {
					forAll(Cg,i) {Cg[i]().storePrevIter();}
					#include "transport/CgiEqn.H"
					}
				}
			
			} // end activate transport
		}
		
		/*print gases 
		if (ph_gcomp>0) { // there are gases only when reaction is present, the freak.g transmit pressures in bars
			Info<<"gvol, Cgs "<<endl;
			for (j=0; j<15;j++) {Info<<gvol[j]<<" "; forAll(Cg,i) {Info<<Cg[i]()[j]<<" ";} Info<<endl;}
		}*/

		//***************  solve reaction  *************************
		// find where the transported conc have changed to calculate only there
		//Info<<"runtime "<<runTime.value()-oldTime<<endl;
		tcnt++;
		if (rSteps<0) {if (tcnt>-rSteps-1) {tcnt=0;} }
		if (rSteps>0) {if (rcnt>rSteps) {rcnt=1;} }

		//if (activateReaction==1 && tcnt==rSteps-1)
		if (rSteps>0) {if (runTime.value() >= wTimes[itwstep-1]+reactStep*rcnt) flg=1;} // here the time to write is a portion of current time period
		if (rSteps<0) {if (tcnt == -rSteps-1) flg=1;} // here the time is a number of flow/transport time steps
		Info<<"for react tcnt "<<tcnt<<" wtime "<<wTimes[itwstep-1]<<" lim "<<wTimes[itwstep-1]+reactStep*rcnt<<" rstep "<<rSteps<<" tcnt "<<tcnt<<" rcnt "<<rcnt<<" flg "<<flg<<endl;
		if (activateReaction==1 && flg)
		{
			rcnt++;flg=0;
			//if (tcnt>rSteps) {tcnt=0;}
			// ***************** find the cells where the conc has changed to calculate there
			//deltaTchem = transportProperties.lookupOrDefault<scalar>("deltaTchem",86400);Info<<"dtchem in reac "<<deltaTchem<<endl;
			std::vector<double> rchange(nxyz,1.);
			if (ncell>1000) {
				for (j=0; j<nxyz;j++) {
					int rj=ractive[j];
					for (i=4; i<ph_ncomp; i++)
					{
						if (abs(c_ph[i*nxyz+j]-Cw[i]()[rj])/(c_ph[i*nxyz+j]+1e-20)<1e-5 || sw[rj]<sw_min[rj]) {rchange[j] = 0.;} else {rchange[j] = 1.;}
					} 
					for (i=0; i<ph_gcomp; i++)
					{
						if (abs(gm_ph[i*nxyz+j]-Cg[i]()[rj])/(gm_ph[i*nxyz+j]+1e-20)<1e-5 || sw[rj]<sw_min[rj]) {rchange[j] = 0.;} else {rchange[j] = 1.;}
					}
				}
			}
			
			// ********************** transfer to phreeqc for c and g, for g we send moles and not pressures, phreeqc does not know the cell volume it considers to be one
			//int icnt = 0;
			forAll(Cw,i) { for (j=0; j<nxyz;j++) {c_ph[i*nxyz+j] = Cw[i]()[ractive[j]];} } 
			forAll(Cg,i) {
					//Vmol[j] = 24.47*(273.15+T()[j])/(273.15+25.);
				for (j=0; j<nxyz;j++) { 
					gm_ph[i*nxyz+j] = max(Cg[i]()[j]*gvol[j],1e-16);// Cg in mol/Lgaz, gm in mol/RV
					} //Cg in moles/L and Vm in L/mol
				} 
			//for (j=0;j<8;j++) {Info<<"p "<<p[j]<<" sw "<<sw[j]<<" gvol "<<gvol[j]<<" T "<<T[j]<<" Cg "; for (i=0; i<ph_gcomp;i++){Info<<Cg[i]()[j]<<" "<<gm_ph[i*nxyz+j]<<" ";} Info<<endl;}
			Info<<" phqVm[0] "<<Vmol[0]<<" "<<endl;
			//auto start = std::chrono::high_resolution_clock::now();
			
			//################# RUN PHREEQC   ################
			// set saturations using rchange
			t_ph.resize(nxyz); p_ph.resize(nxyz);// showul dnot be necessary but seems required???
			Info<<"rchange ";
			for (j=0; j<nxyz;j++) {rchange[j] *= sw[j];t_ph[j]=T()[j];if(j<10) {Info<<" "<<rchange[j];} } Info<<endl;
			//for (j=0; j<nxyz;j++) {p_ph[j]=p[j]/atmPa;}
			if (ph_gcomp>0) {freak.setGvol(gvol);} // set gas volume in phreeqc
			freak.setWsat(rchange); // rchange for the calculation domain, with 0 outside, sw saturation
			freak.setC(c_ph);//transfer c_ph to freak : it does not work to send directly to freak.c
			freak.setGm(gm_ph);//transfer gm_ph to freak
			if (activateThermal==1) {freak.setTemp(t_ph);}
			//freak.setP(p_ph);//transfer pressure to freak
			freak.setTstep((runTime.value()-oldTimeReac)*tunits); //Info<<" this tme "<< runTime.value()<<" old "<<oldTime<<endl;//the calculation time shall include all time since las phreeqc run
			Info << "running phreeqc dt "<<(runTime.value()-oldTimeReac)*tunits<<endl;
			int a0= phqRun(freak);Info<<"end results "<<a0<<endl;
			if (a0==-7) { // irm _fail : what to do
				forAll(Cw,ic) { if (immobile[ic]==0) {Cw[ic]()=Cw[ic]().prevIter();} } // back to previous conc values
				tcnt=rSteps-2;//insure that phreeqc will be run on next time step !!!!! oldTimeReac to be corrected
				rewind=1;Info<<endl;continue;
				} 
			a0=getSelOutput(freak);
			Info << "phreeqc done "<<endl;

			//auto finish = std::chrono::high_resolution_clock::now();
			//std::chrono::duration<double> dt = finish - start;dure = dure+dt.count();
			
			// #############  transfer back to C 
			
			forAll(Cw,ic) {
				// if transfered as source term for mobile species, but immobile stay here
				if (immobile[ic]==1) {
					for (j=0;j<nxyz;j++) {Cw[ic]()[ractive[j]] = freak.c[ic*nxyz+j];} // transfer phq -> opf		
				}
				
				//for (j=0;j<nxyz;j++) {if (bcCwi[ractive[j]]==0) {Cw[ic]()[ractive[j]] = freak.c[ic*nxyz+j];} }// transfer phq -> opf		
			}
			
			//Info<<"Cw mid "<<Cw[7]()[int(ncell/2)]<<endl;

			//----------- finding the time step (not for immobile species)
			if (reactionSteps==0) // in this case, we use dtForChem
				{
				dtForChem = 1e12;
				forAll(Cw,ic) // dissolved
					if ((ic>3)&&(immobile[ic]==0))
					{
						dC = 0;dff=0;mxC=0;mnC=1;
						for (j=0;j<nxyz;j++) 
							{mnC=min(mnC,freak.c[ic*nxyz+j]);
							mxC=max(mxC,freak.c[ic*nxyz+j]);}
						for (j=0; j<nxyz;j++)
							if (bcCwi[j]==0) 
							{
							dff = mag(Cw[ic]()[ractive[j]]-freak.c[ic*nxyz+j]);
							dC = max(dC,dff);
							//Info<<"ic, j "<<ic<<" "<<j<<" c "<<Cw[ic]()[ractive[j]]<<" "<<freak.c[ic*nxyz+j]<<" "<<dff<<" "<<dC<<endl;
							}
						dC = dC/(mxC-mnC+SMALL);
						dtForChem = min(dtForChem,dCmax/(max(dC,0)+SMALL)*runTime.deltaTValue()); Info<<"ic "<<ic<<" dC "<<dC<<" dtForChem "<<dtForChem<<endl;
					} 
				Info<< "dt "<<runTime.deltaTValue()<<" dC "<<dC<<" dtForC "<<dtForC<<" dtForChem " << dtForChem << endl; 
				newDeltaT = min(dtForChem, newDeltaT);
				if ((dtForChem>dtForC*10)&&(tstep>10)) {rSteps=10;} //even if dtForChem is high, keep however a mini value of for security
				else if ((dtForChem>dtForC*2)&&(tstep>10)) {rSteps=static_cast<int>(std::round(dtForChem/dtForC));}
				else {rSteps=1;}
				}
			Info<<"new deltat "<<newDeltaT<<endl;

			//************************ gas, read partial pressures (freak.g in atm) set Cg to fraction (freak.g/Cgtot) and set p to sum of Cg
			
			if (ph_gcomp>0) // multiphase
				{
				float r;
				for (j=0; j<nxyz;j++) //should consider ractive
					{
					//Cgtot = 0; Gmtot = 0;
					//forAll(Cg,i) {freak.g[i*nxyz+j] = max(freak.g[i*nxyz+j],0.); Cgtot += freak.g[i*nxyz+j]; Gmtot += freak.gm[i*nxyz+j];}
					forAll(Cg,i) {Cg[i]()[j] = freak.gm[i*nxyz+j]/gvol[j];} 
					if ((flowType==4)&&(activateEbullition==1)){
						//Vmol[j] = freak.gvol[j]/Gmtot; //1.003 strange but needed to have equilibrated pressure
						//r = freak.p[j]/(p[j]/atmPa);//change of gas pressure
						//sw[j]=1-r*freak.gvol[j]/eps[j];
						p[j] = freak.p[j]*atmPa; //.g are gas pressures, tot is the sum
						}
					//sw[j] = freak.wsat[j];
					//if (j<5) {Info<<"p "<<p[j]<<" sw "<<sw[j]<<" gvol "<<gvol[j]<<" Cg "; for (i=0; i<ph_gcomp;i++){Info<<Cg[i]()[j]<<" ";} Info<<endl;}
					}
				}
			
			// find the variation of wsat from nb moles H2O in gas phase, phreeqc considers a volume of 1 dm3
			if (freak.iGwater>-1) //should consider ractive
				{
				double a1=0.;
				//iw = freak.iGwater; done at start
				for (j=0; j<nxyz;j++)
					{
					a1 = max(freak.gm[iw*nxyz+j],0.) - gm_ph[iw*nxyz+j]; // delta gm water in moles in gvol (equiv of 1L of medium volume)
					//a1 = a1 *   //gm_ph = Cg[i]()*gvol/Vmol
					//if (j<5) {Info<< " gm_ph "<< gm_ph[iw*nxyz+j] << " frk "<< freak.gm[iw*nxyz+j] <<" a1 "<<a1<< endl;}
					sw[j] = max(sw_min[j],sw[j] - a1*.01801/eps[j]); 
					}
				} 
				//nb of moles of H2O(g) transformed in water volume (1 mol 18.01 mL at 25°C)
			//for (j=0;j<3;j++) {Info <<" new sw "<<sw[j]<< endl;}
				
			oldTimeReac = runTime.value()*1;
			reactIndicator = 1; //to know after reaction
		} //end activate reaction

		//runTime.setDeltaT(min(newDeltaT,maxDeltaT)); //deltaT must be set only once!!!
		//bool ts = runTime.write();
		//write
		#include "observation.H"
		#include "budget.H"
		
		if (flagW==1) {runTime.writeNow();tcnt=0;Info<<"l548, writing"<<endl;}
		
		//if (flowType==4) {phiGr.write();}
		if (activateReaction==1  && flagW==1) {
			phiw.write();phig.write();
			std::ofstream outFile(cur_dir/runTime.timeName()/"Species");
			outFile.unsetf(std::ios::scientific);outFile.precision(6);
			for (j=0;j<nxyz;j++)
				{ for (i=0;i<freak.nselect;i++) {outFile << freak.spc[i*nxyz + j]<<" ";} outFile <<"\n"; }
			}
			
		if (flagW==1) {flagW=0;}
		Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
			<< "  ClockTime = " << runTime.elapsedClockTime() << " s"
			<< nl << endl;

		Info<< "End tstep\n" << endl;

	}
	
	runTime.writeNow();
	/*
	if (activateReaction==1) {
		phiw.write();phig.write();
		std::ofstream outFile(cur_dir/runTime.timeName()/"Species");
		outFile.unsetf(std::ios::scientific);outFile.precision(6);
		for (j=0;j<nxyz;j++)
			{ for (i=0;i<freak.nselect;i++) {outFile << freak.spc[i*nxyz + j]<<" ";} outFile <<"\n"; }
		}
	*/
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
