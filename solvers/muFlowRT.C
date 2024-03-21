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
    darcyFoam

Description
    Stationary solver for incompressible single-phase flow in porous medium

Developers
    Pierre Horgue then Olivier Atteia to integrate phreeqc

\*---------------------------------------------------------------------------*/
//#include <cmath.h>
//#include <iomanip> //NB when in < > don't add the .h
#include <stdlib.h>
#include <vector>
//#include <iostream>
#include <fstream>
#include <iterator>
#include <sstream>
#include <string>

#include <chrono>  // for high_resolution_clock
#include "fvCFD.H"
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
std::vector<double> c_ph,gm_ph,t_ph,p_ph,gvol,ractive,solu_conc,gas_conc;
std::vector<int> immobile,wTimes;
float atmPa=101325.;float vmw,Cgtot,Gmtot;
int i,j,iw,oindex;
my_phq freak; //need to be here to be availabel for every chem condition

#include "myFunc.H"
// read the myfunc file
dicFunc fDe_T;		   		  

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;
#include "utilities.h" // for reading binary reading tables..
#include "plugin_Cgi.H"
plugin_Cgi plugCgi;


int main(int argc, char *argv[])
{
	//init openFoam
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
	//#include "polyMesh.H" was supposed to be usefull for cellLevles but creates pb elsewhere
    #include "readGravitationalAcceleration.H"
    #include "createFields.H"
	#include "readPicardControls.H"
    simpleControl simple(mesh);
	#include "createThetaFields.H"
	#include "create2phaseFields.H"
	//parms of the De=f(T) function
	#include "readFunc.H"
	Info << "fDe parms "<< fDe_T.fparms[0] << " "<<fDe_T.fparms[1] << endl;
	#include "transport/createCFields.H"
	#include "transport/createTFields.H"
	
	// reading some general files (times, obs points)
	std::ifstream inputwTimes{cur_dir+"/constant/options/writetimes" }; // version 0 shall contain 0 for inactive and 1 for active reaction cell
	wTimes = {std::istream_iterator<int>{inputwTimes}, {}};
	
	//Observations : a file with the name of each pooint and its x,y,z, coordinates
	fname = cur_dir+"/constant/options/obspts";std::vector<float> obs;int nobs;std::vector<int> icello;outTable observ;
	if (fexists(fname))
		{
		observ = readTable(fname);nobs = observ.nrow;int ncol = observ.ncol;icello.resize(nobs);std::cout<<"nobs "<<nobs<<" "<<ncol<<"\n";
		for (int io=0;io<nobs;io++) {
			//vector coord(observ.data[io*ncol],observ.data[io*ncol+1],observ.data[io*ncol+2]);
			int ix=observ.data[io*ncol];int iz=observ.data[io*ncol+1];
			icello[io]=iz*ncell_lay+ix;std::cout<<"obs i "<<ix<<" "<<iz<<" "<<icello[io]<<"\n";
			}
		}		
		
	if (activateReaction==1)
	{
	//##############  phreeqc intiialisation for solutions and gases
	std::ifstream inputRactive{cur_dir+"/constant/options/ractive" }; // version 0 shall contain 0 for inactive and 1 for active reaction cell
	ractive = {std::istream_iterator<int>{inputRactive}, {}};
	std::ifstream inputInit{cur_dir/"phqinit.txt"};
	std::vector<int> ph_init{std::istream_iterator<int>{inputInit}, {}}; //for (int i=0; i<7;i++){Info << "init nb "<< ph_data[i] << endl;}
	nxyz=ph_init[0];ph_ncomp=ph_init[1];ph_gcomp=ph_init[2];ph_nsolu=ph_init[3];
	freak.setDB(cur_dir/"phreeqc.dat");
	//make a first init to calculate the inital solutions (poro required for the initial solutions)
	freak.setChemFile(cur_dir/"initChem.pqi"); //Info << "initCh read " << endl;
	freak.setData(ph_init);
	std::vector<double> poro(ph_nsolu,0.25);std::vector<double> t_ph(ph_nsolu,0.);
	for (i=0;i<ph_nsolu;i++) {poro[i] = eps[i];t_ph[i]=T[i];}
	//gvol.resize(ph_nsolu,1.); 
	//for (i=0;i<ph_nsolu;i++) {gvol[i] = eps[i]*0.1;}
	//p_ph.resize(ph_nsolu,1.);
	//std::vector<double> wsat(ph_nsolu,0.9); // at teh beginning we set high gaz volume so the solution does not modify the gaz composition
	freak.setPoro(poro);freak.setTemp(t_ph);
	freak.init(); //Info << "nxyz " << nxyz << endl;
	
	//################ writes the initial solutions and gases to files
	std::ofstream outFile(cur_dir/"constant/options/solutions");
	solu_conc.resize(ph_nsolu*ph_ncomp,0.);Info << "nsolu "<<ph_nsolu << " ncomp "<< ph_ncomp <<endl;
	gas_conc.resize(ph_nsolu*ph_gcomp,0.);
	for (i=0;i<ph_nsolu;i++) // solu number
		{ for (j=0;j<ph_ncomp;j++) // component number
			{
			float a = freak.c[j*ph_nsolu+i];
			solu_conc[i*ph_ncomp+j] = a; outFile << a << "\n"; 
			} 
		}
	outFile.close();
	
	std::ofstream outFile1(cur_dir/"constant/options/gases");
	for (i=0;i<ph_nsolu;i++) // solu number
		{ 
		Cgtot = 0;
		for (j=0;j<ph_gcomp;j++) // loop component to calculate Cgtot
			{ Cgtot += freak.g[j*ph_nsolu+i];Info<<"Cgtot "<<Cgtot<<endl;}
		for (j=0;j<ph_gcomp;j++) // loop over component to write Cgi		
			{ float a = freak.g[j*ph_nsolu+i]/Cgtot; 
			gas_conc[i*ph_gcomp+j] = a;outFile1 << a << "\n"; } //in fraction /Cgtot/phreeqcVm
		}
	outFile1.close();

	//##################### make the initialization for the full domain in phreeqc : data,poro, gvol
	std::ifstream inputData{cur_dir/"phqfoam.txt"};
	std::vector<int> ph_data{std::istream_iterator<int>{inputData}, {}}; //for (int i=0; i<7;i++){Info << "init nb "<< ph_data[i] << endl;}
	nxyz=ph_data[0];ph_ncomp=ph_data[1];ph_gcomp=ph_data[2];ph_nsolu=ph_data[3]; //!!! nxyz here is inside the ractive part
	freak.setData(ph_data); Info << "nxyz " << nxyz << endl;
	//initiate poro and gas volume
	poro.resize(nxyz,0); t_ph.resize(nxyz,0.);
	for (i=0;i<nxyz;i++) {poro[i]=eps[i];t_ph[i]=T[i];}
	//wsat.resize(nxyz,0.994);
	freak.setPoro(poro);freak.setTemp(t_ph);
	//freak.setWsat(wsat); // rchange for the calculation doamin, with 0 outside, sw saturation
	p_ph.resize(nxyz);
	for (i=0;i<nxyz;i++) {p_ph[i]=p[i]/atmPa;}	
	freak.init();
	
	//##############" build the c_ph and gm_ph fields and get conc from phreeqc (c_ph=freak.c but needed two variables for format questions)
	c_ph.resize(nxyz*ph_ncomp,0);
	for (i=0; i<ph_ncomp;i++)
		for (int j=0;j<nxyz;j++) 
			{c_ph[i*nxyz+j] = freak.c[i*nxyz+j];}
			// gases are in atm (freak.g) we start with ideal gas
	gm_ph.resize(nxyz*ph_gcomp,1e-9);gvol.resize(nxyz);phreeqcVm.resize(nxyz); // moles of gas
	//in order to have correct correct pressure, we set gm using only the frist gas
	if (ph_gcomp>0)
		{
		for (int j=0;j<nxyz;j++) 
			{
			phreeqcVm[j] = 24.47*(273.15+25.)/(273.15+25.)/(p[j]/atmPa);
			//gm_ph[i*nxyz+j] = freak.gm[i*nxyz+j];
			gvol[j] = max(eps[j]*(1-sw[j]),1e-4);
			gm_ph[j] = max(gvol[j]/phreeqcVm[j],1e-16); //only first gas comp
			}
		freak.setGvol(gvol); // set gas volume in phreeqc
		freak.setGm(gm_ph);//transfer gm_ph to freak
		freak.setTemp(t_ph);
		freak.setP(p_ph);
		freak.run();
		for (int j=0;j<nxyz;j++)  {std::cout<<"p t0 "<<freak.p[j]<<"gv t0 "<<freak.gvol[j]<<" Vm t0 "<<phreeqcVm[j]<<"\n";}
		}
	
	iw = freak.iGwater;	
	

	Info<<"end c_ph and gm_ph "<<endl;
	//###################"" loading immobile component is present  ###############""
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
		ph_ncomp=0;ph_gcomp=0;
	}
	
	plugCgi.init(cur_dir,transportProperties,mesh,freak); // initiate the plugin for transport

	// #include "createFvOptions.H"
	#include "transport/createCwiFields.H"
	#include "transport/createCgiFields.H"

	//################### attribute to Cwi and Cgi concentrations/pressures from phreeqc
	int icnt = 0;
	//Info<<"ractive 2 "<<ractive[2]<<endl;
	//Info<<" gmph 2 "<<  gm_ph[2]<<" frk.c size "<<freak.c.size()<<endl;
	
	if (activateReaction==1) 
	{
		forAll(Cw,i) {
			for (j=0; j<nxyz;j++){Cw[i]()[ractive[j]] = freak.c[i*nxyz+j];} // transfer freak.c to Cw
			} 
		if (ph_gcomp>0) { // there are gases only when reaction is present, the freak.g transmit pressures in bars
			for (j=0; j<nxyz;j++)
				{
				Cgtot = 0;
				forAll(Cg,i){Cgtot += freak.g[i*nxyz+j];} 
				forAll(Cg,i){Cg[i]()[j] = freak.g[i*nxyz+j]/Cgtot;} // /Cgtot/phreeqcVm;}
				}
		}
	}
	if (ph_gcomp>1) {Info<<" cg 0 1 "<<Cg[0]()[1]<<endl;}
		
	
	//######################## run the steady state for hp
	dimensionedScalar st = runTime.startTime();
	dimensionedScalar et = runTime.endTime();
	dimensionedScalar dt0 = mesh.time().deltaTValue();
	scalar dt1 = runTime.controlDict().lookupOrDefault("writeInterval",0)/10;Info<<"dt1 "<<dt1<<endl;
	scalar residu;
	if ((flowStartSteady==1)&&(flowType>0)&&(flowType<=2))
		{
		runTime.setDeltaT(1);
		while (simple.loop(runTime))
			{
			#include "hstdEqn.H"
			}
		}
	Info <<"st time "<<st<<endl;
	//runTime.setTime(st,0); // 12/3/21 time value and index
	runTime.setEndTime(et); Info<<"end "<<runTime.endTime()<<endl;
	runTime.runTimeModifiable();
	//runTime.read();
	Info <<"dt time "<<dt0<<endl;
	runTime.setDeltaT(dt0);
	float oldTime=0;
	Info<<"time rebuilt st "<<runTime.startTime()<<" dt "<<runTime.deltaTValue()<<endl;
	
	//const IOobject meshObj("constant/polyMesh", runTime.constant(),"volScalarField", "h", IOobject::MUST_READ);
    //const fvMesh mesh(meshObj);

    // Récupération des cellules réordonnées
    //const labelList& cellLevel = mesh.cellLevels();

// Obtention des niveaux de raffinement pour la cellule 42
//const labelList& localIndices = mesh.cells().local();
    // Affichage des index des cellules réordonnées
    //Info << "Index des cellules réordonnées : " << level << endl;
	
	int itstep = 0;int tcnt = 0;float wtime=wTimes[itstep];int flagW=0; // first tstep is 0
	while (runTime.run())
    {
		//double rt = runTime.controlDict().lookupOrDefault("writeInterval",0);
		//set time step to stop at writeTimes
		oldTime = mesh.time().value();
		float dt=wtime - oldTime;
		if ((dt <= float(runTime.deltaTValue())*(1+5e-5))||(std::round(dt)==std::round(runTime.deltaTValue()))) //wpb round jus tone, sometime for long times there is a diff of 2 or 3 sec
			{
			Info<<"dt "<<dt<<" deltaT "<<runTime.deltaTValue();
			//float deltaTFact = dt/float(runTime.deltaTValue())*(1+1e-3);
			runTime.setDeltaTNoAdjust(dt);itstep+=1;wtime=wTimes[itstep];flagW=1;
			Info<<" flg "<<flagW<<endl;}
		Info<<"i time "<<itstep<<" oldt "<<oldTime<<" wt "<<wtime<<" dt "<<dt<<" deltaT "<<float(runTime.deltaTValue())+0.01<<" flg "<<flagW<<endl;
		if (dt==0) {itstep+=1;wtime=wTimes[itstep];flagW=1;}
		Info<<"i time "<<itstep<<" "<<wtime<<" "<<flagW<<endl;
		runTime.read();
		// #include "transport/setDeltaTtrsp.H"
		runTime++;
		//if (mesh.time().value()==wtime) {flagW=1;}
		Info << "time = " << runTime.timeName() <<  "  deltaT = " <<  runTime.deltaTValue() << endl;
		// *********** here provide change of density and viscosity if required
		
		//***********************  solve transient flow   *******************************
		//Info<<"p before flow ";for (int j=0; j<nxyz;j++) {Info<<p[j]/atmPa<<" ";};Info<<endl;
		//Info<<"sw before flow ";for (int j=0; j<nxyz;j++) {Info<<sw[j]<<" ";};Info<<endl;
		if (flowType>0) {
			#include "flow.H"
			}
		//Info<<"p after flow ";for (int j=0; j<nxyz;j++) {Info<<p[j]/atmPa<<" ";};Info<<endl;
		//Info<<"sw after flow ";for (int j=0; j<nxyz;j++) {Info<<sw[j]<<" ";};Info<<endl;
		if (ph_gcomp>0) {for (j=0; j<nxyz;j++) {gvol[j]=max(eps[j]*(1-sw[j]),1e-4);} }
		
		//***************  solve Transport  *************************
		if (activateThermal==1) {
			#include "transport/TEqn.H"
			}
		if (activateTransport==1) {
			if (activateReaction==0) {
				//#include "transport/setDeltaTtrsp.H"
				#include "transport/CEqn.H"
				}
			
			else {
				forAll(Cw,i) {Cw[i]().storePrevIter();} // for cells outside calculation
				#include "transport/CwiEqn.H"
				if (ph_gcomp>0) {
					forAll(Cg,i) {Cg[i]().storePrevIter();}
					#include "transport/CgiEqn.H"
					}
				}
			
			}
		runTime.setDeltaT (min(newDeltaT,maxDeltaT)); //deltaT msut be set only once!!!

		//***************  solve reaction  *************************
		// find where the transported conc have changed to calculate only there
		//Info<<"runtime "<<runTime.value()-oldTime<<endl;
		tcnt++;
		if (tcnt>reactionSteps-1) {tcnt=0;}
		//if (activateReaction==1 && deltaTchem<=0)  //runTime.value()-oldTime>dT1
		
		if (activateReaction==1 && tcnt==reactionSteps-1)
		{
			// find the cells where the chcemistry has changed to calculate there
			//icnt = 0;
			deltaTchem = transportProperties.lookupOrDefault<scalar>("deltaTchem",86400);Info<<"dtchem in reac "<<deltaTchem<<endl;
			std::vector<double> rchange(nxyz,0.);
			for (i=4; i<ph_ncomp; i++)
				{
				for (j=0; j<nxyz;j++)
					{
					if (abs(c_ph[i*nxyz+j]-Cw[i]()[ractive[j]])/(c_ph[i*nxyz+j]+1e-20)>1e-4 && sw[ractive[j]]>sw_min[j]*1.5) {icnt++;rchange[j] = 1.;}
					}
				} 
			for (i=0; i<ph_gcomp; i++)
				{
				for (j=0; j<nxyz;j++)
					{
					if (abs(gm_ph[i*nxyz+j]-Cg[i]()[ractive[j]])/(gm_ph[i*nxyz+j]+1e-20)>1e-5 && sw[ractive[j]]>sw_min[j]) {rchange[j] = max(rchange[j],1.);}
					}
				} 
			// transfer to phreeqc for c and g, for g we send moles and not pressures, phreeqc does not know the cell volume it considers to be one
			//int icnt = 0;
			for (j=0; j<nxyz;j++) { forAll(Cw,i) {c_ph[i*nxyz+j] = Cw[i]()[ractive[j]];} } 
			for (j=0; j<nxyz;j++) { 
				forAll(Cg,i) {
					if (flowType!=4) {phreeqcVm[j] = 24.47*(273.15+T()[j])/(273.15+25.);}
					gm_ph[i*nxyz+j] = max(Cg[i]()[j]*gvol[j]/phreeqcVm[j],1e-16);
					} //Cg in fraction and Vm in L/mol
				} 

			
			Info<<" phqVm[0] "<<phreeqcVm[0]<<" "<<endl;
			
			//auto start = std::chrono::high_resolution_clock::now();
			//################# RUN PHREEQC   ################
			// set saturations using rchange
			t_ph.resize(nxyz); p_ph.resize(nxyz);// showul dnot be necessary but seems required???
			for (j=0; j<nxyz;j++) {rchange[j] *= sw[j];t_ph[j]=T()[j];} //Info<<" "<<rchange[j];
			for (j=0; j<nxyz;j++) {p_ph[j]=p[j]/atmPa;}
			freak.setGvol(gvol); // set gas volume in phreeqc
			freak.setWsat(rchange); // rchange for the calculation doamin, with 0 outside, sw saturation
			freak.setC(c_ph);//transfer c_ph to freak : it does not work to send directly to freak.c
			freak.setGm(gm_ph);//transfer gm_ph to freak
			freak.setTemp(t_ph);
			freak.setP(p_ph);//transfer pressure to freak
			freak.setTstep(runTime.value()-oldTime); //Info<<" this tme "<< runTime.value()<<" old "<<oldTime<<endl;//the calculation time shall include all time since las phreeqc run
			Info << "running phreeqc dt "<<runTime.value()-oldTime<<endl;
			freak.run();
			freak.getSelOutput();
			Info << "phreeqc done "<<endl;
				
			//auto finish = std::chrono::high_resolution_clock::now();
			//std::chrono::duration<double> dt = finish - start;dure = dure+dt.count();
			
			// transfer back to C but before calc tstep
			
			// set to previous values (outisde rchange)
			forAll(Cw,i) {Cw[i]() = Cw[i]().prevIter();}
			if (ph_gcomp>0) {forAll(Cg,i) {Cg[i]() = Cg[i]().prevIter();} }

			
			forAll(Cw,ic) // dissolved
				{
					dC0 = 0;
					for (j=0; j<nxyz;j++)
						{
						Cw[ic]()[ractive[j]] = freak.c[ic*nxyz+j];
						if (j==4) {Info<<"ic "<<i<<" c "<<Cw[ic]()[4]<<endl;}
						}
				} 
			// gas, read partial pressures (freak.g in atm) set Cg to fraction (freak.g/Cgtot) and set p to sum of Cg
			if (ph_gcomp>0) // multiphase
				{
				for (j=0; j<nxyz;j++) //should consider ractive
					{
					Cgtot = 0; Gmtot = 0;
					forAll(Cg,i) {freak.g[i*nxyz+j] = max(freak.g[i*nxyz+j],0.); Cgtot += freak.g[i*nxyz+j]; Gmtot += freak.gm[i*nxyz+j];}
					forAll(Cg,i) {Cg[i]()[j] = freak.g[i*nxyz+j]/Cgtot;}
					if (flowType==4) {
						gvol[j]=Gmtot*phreeqcVm[j];phreeqcVm[j] = gvol[j]/Gmtot;
						//if (gvol[j]>1.01e-4) 
						sw[j]=1-gvol[j]/eps[j];
						}//p[j] = Cgtot*atmPa;}
					//sw[j] = freak.wsat[j];
					}
				for (i=0; i<ph_gcomp;i++){Info <<"g_spc 0 "<< i <<" Cg "<< Cg[i]()[0] <<" p "<<p[0]<<" sw "<<sw[0]<< endl;}
				for (i=0; i<ph_gcomp;i++){Info <<"g_spc 20 "<< i <<" Cg "<< Cg[i]()[20] <<" p "<<p[20]<<" sw "<<sw[20]<< endl;}
				for (i=0; i<ph_gcomp;i++){Info <<"g_spc 50 "<< i <<" Cg "<< Cg[i]()[50] <<" p "<<p[50]<<" sw "<<sw[50]<< endl;}
				}
			// find the variation of wsat from nb moles H2O in gas phase, phreeqc considers a volume of 1 dm3
			if (freak.iGwater>-1) //should consider ractive
				{
				double a1=0.;
				//iw = freak.iGwater; done at start
				for (j=0; j<nxyz;j++)
					{
					a1 = max(freak.gm[iw*nxyz+j],0.) - gm_ph[iw*nxyz+j]; // delta gm water in moles in gvol (equiv of 1L of medium volume)
					//a1 = a1 *   //gm_ph = Cg[i]()*gvol/phreeqcVm
					if (j<5) {Info<< " gm_ph "<< gm_ph[iw*nxyz+j] << " frk "<< freak.gm[iw*nxyz+j] <<" a1 "<<a1<< endl;}
					sw[j] = max(sw_min[j],sw[j] - a1*.01801/eps[j]); 
					}
				} 
				//nb of moles of H2O(g) transformed in water volume (1 mol 18.01 mL at 25°C)
			for (j=0;j<3;j++) {Info <<" new sw "<<sw[j]<< endl;}
				
			oldTime = runTime.value()*1;
		} //end activate reaction
		
		
		//bool ts = runTime.write();
		//write
		#include "observation.H"
		#include "budget.H"
		if (flagW==1) {runTime.writeNow();}
		
		//if (flowType==4) {phiGr.write();}
		if (activateReaction==1  && flagW==1) {
			phiw.write();phig.write();
			std::ofstream outFile(cur_dir/runTime.timeName()/"Species");
			outFile.unsetf(std::ios::scientific);outFile.precision(6);
			for (j=0;j<nxyz;j++)
				{
				for (i=0;i<freak.nselect;i++) {outFile << freak.spc[i*nxyz + j]<<" ";}
				outFile <<"\n";
				}
			}
			
		if (flagW==1) {flagW=0;}
		Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
			<< "  ClockTime = " << runTime.elapsedClockTime() << " s"
			<< nl << endl;

		Info<< "End tstep\n" << endl;

	}
    return 0;
}

// ************************************************************************* //
