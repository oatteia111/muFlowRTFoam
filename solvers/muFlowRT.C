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
#include <iostream>
#include <fstream>
#include <iterator>
#include <sstream>
#include <string>

#include <chrono>  // for high_resolution_clock
#include "fvCFD.H"
#include "incompressiblePhase.H"
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


#include "phreeqc/initPhreeqc.H"

std::vector<double> a(12,0.);
std::vector<double> c_ph,gm_ph,t_ph,p_ph,gvol,ractive,solu_conc,gas_conc;
float atmPa=101325.;float vmw,Cgtot,Gmtot;
int i,j,iw;
my_phq freak; //need to be here to be availabel for every chem condition
#include "myFunc.H"
// read the myfunc file
dicFunc fDe_T;		   		  

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;

std::vector<int> indexC(labelList &cells, std::vector<float> &data)
{
    std::vector<int> c1(cells.size(),0);//Info<<"in index "<<cells.size()<<endl;
	for (int i=0; i<cells.size();i++)  // reads the first ncells lines
		{
		auto iter = std::find(cells.begin(), cells.end(), static_cast<int>(data[i*4+1]));
		int i1 = {std::distance(cells.begin(), iter)};  //Info << i << " icd "<< icd << " cll " << cells_[i] << " indx " << a << endl; // cell number in cellsData //index of cellsData in cells_
		c1[i] = cells[i1];
		}
    return c1;
}

inline bool fexists(const std::string& name) {
    ifstream f(name.c_str());
    return f.good();
}

struct outData {float t; std::vector<float> d;};
 
// a function to get data from binary file
outData getCbuffer(char* fname, int itime, int ncell) {
    std::vector<float> data(ncell*4);
	std::ifstream inputData{cur_dir+"/constant/options/"+fname, std::ios::binary}; //
	//while (inputHdrndata.read(reinterpret_cast<char*>(&f0), sizeof(float))) {cellsHdrnData.push_back(f0);}
	//for (int ic=0;ic<ncell;ic++) {inputData.read(reinterpret_cast<char*>(&f0), sizeof(float));data.push_back(f0);}
	// position at the place
	inputData.seekg(ncell*itime*4*sizeof(float)); //each line is composed of 4 numbers and there are two values at the beginning

    inputData.read(reinterpret_cast<char*>(&data[0]), ncell*4*sizeof(float));
	float time;	
	inputData.read(reinterpret_cast<char*>(&time), sizeof(float));
	// other possiblity
	/*char buffer[BUFFERSIZE];
	FILE * filp = fopen("filename.bin", "rb"); 
	int bytes_read = fread(buffer, sizeof(char), BUFFERSIZE, filp); */
	outData output;
	output.t = time;std::cout<<"readbin "<<fname <<" "<<itime<<" "<<ncell;
	output.d = data;
    return output;
}

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
	std::vector<double> parm0(1,1e-3); // the parameters for the modification of DEffg with temperature
	#include "readFunc.H"
	Info << "fDe parms "<< fDe_T.parms[1] << endl;																																   
	
	if (activateReaction==1)
	{
	//##############  phreeqc intiialisation for solutions and gases
	float rv=1.; // dont change it, it does not work
	std::ifstream inputRactive{cur_dir+"/constant/options/ractive" }; // version 0 shall contain 0 for inactive and 1 for active reaction cell
	ractive = {std::istream_iterator<int>{inputRactive}, {}};
	std::ifstream inputInit{cur_dir/"phqinit.txt"};
	std::vector<int> ph_init{std::istream_iterator<int>{inputInit}, {}}; //for (int i=0; i<7;i++){Info << "init nb "<< ph_data[i] << endl;}
	nxyz=ph_init[0];ph_ncomp=ph_init[1];ph_gcomp=ph_init[2];ph_nsolu=ph_init[3];
	freak.setDB(cur_dir/"phreeqc.dat");
	//make a first init to calculate the inital solutions (poro required for the initial solutions)
	freak.setChemFile(cur_dir/"initChem.pqi"); //Info << "initCh read " << endl;
	freak.setData(ph_init);
	std::vector<double> poro(ph_nsolu,0.25);
	for (i=0;i<ph_nsolu;i++) {poro[i] = eps[i];}
	//gvol.resize(ph_nsolu,1.); 
	//for (i=0;i<ph_nsolu;i++) {gvol[i] = eps[i]*0.1;}
	//p_ph.resize(ph_nsolu,1.);
	//std::vector<double> wsat(ph_nsolu,0.9); // at teh beginning we set high gaz volume so the solution does not modify the gaz composition
	freak.setPoro(poro);
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
	poro.resize(nxyz,0);
	for (i=0;i<nxyz;i++) {poro[i]=eps[i];}
	//wsat.resize(nxyz,0.994);
	freak.setPoro(poro);
	//freak.setWsat(wsat); // rchange for the calculation doamin, with 0 outside, sw saturation
	//for (i=0;i<nxyz;i++) {p_ph[i]=p[i]/atmPa;}	
	//freak.setP(p_ph);
	freak.init();
	gvol.resize(nxyz,0.01);
	t_ph.resize(nxyz,20.);
	//freak.run();
	
	//##############" build the c_ph and gm_ph fields and get conc from phreeqc (c_ph=freak.c but needed two variables for format questions)
	c_ph.resize(nxyz*ph_ncomp,0);
	for (i=0; i<ph_ncomp;i++)
		for (int j=0;j<nxyz;j++) 
			{c_ph[i*nxyz+j] = freak.c[i*nxyz+j];}
	// gases are in atm (freak.g) we start with ideal gas
	gm_ph.resize(nxyz*ph_gcomp,0); // moles of gas
	Info<<"gcomp "<<ph_gcomp<<" g size "<<freak.g.size()<<" frk.c size "<<freak.c.size()<<endl;
	for (i=0; i<ph_gcomp;i++)
		for (int j=0;j<nxyz;j++) 
			{
			phreeqcVm[j] = 24.5/(p[j]/atmPa);
			gm_ph[i*nxyz+j] = freak.gm[i*nxyz+j];
			}
	iw = freak.iGwater;
	
	} //end of activateReation
	else  //only flow or flow+transport
	{
		ph_ncomp=0;
	}
	// #include "createFvOptions.H"
	#include "transport/createCFields.H"
	#include "transport/createTFields.H"
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
	dimensionedScalar dt = mesh.time().deltaTValue();
	scalar dt1 = runTime.controlDict().lookupOrDefault("writeInterval",0)/10;Info<<"dt1 "<<dt1<<endl;

	if ((flowStartSteady==1)&&(flowType>0))
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
	Info <<"dt time "<<dt<<endl;
	runTime.setDeltaT(dt);
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
	
	int istep = 0;int tcnt = 0;
	while (runTime.run())
    {
		runTime.read();
		// #include "transport/setDeltaTtrsp.H"
		runTime++;
		Info << "time = " << runTime.timeName() <<  "  deltaT = " <<  runTime.deltaTValue() << endl;
		// *********** here provide change of density and viscosity if required
		
		//***********************  solve transient flow   *******************************
		//for (j=0; j<nxyz;j++) {if (j<6) {Info<<"p before flow "<<p[j]/atmPa<<" sw "<<sw[j]<<endl;}}
		if (flowType>0) {
			#include "flow.H"
			}
		deltaTchem -= runTime.deltaTValue(); Info<<"dtchem "<<deltaTchem<<endl;
		if (ph_gcomp>0) {for (j=0; j<nxyz;j++) {gvol[j]=eps[j]*(1-sw[j]);} }
		
		//***************  solve Transport  *************************
		if (activateThermal==1) {
			#include "transport/TEqn.H"
			}
		if (activateTransport==1) {
			if (activateReaction==0) {
				#include "transport/setDeltaTtrsp.H"
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
		//dC1 = dC*.999;dT1 = dT*.999;Info<<"dC1  "<<dC1<<endl;
		//Info<<" cg 0 0 "<<Cg[0]()[0]<<" cg 0 1 "<<Cg[0]()[1]<<endl;
		//if (ph_gcomp>1) {Info<<" cg 0 1 "<<Cg[0]()[1]<<endl;}
		
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
			if (flowType == 4)
				{
				for (j=0; j<nxyz;j++) { 
					forAll(Cg,i) {
					gm_ph[i*nxyz+j] = max(Cg[i]()[j]*gvol[j]/phreeqcVm[j],1e-6);//Cg in fraction and Vm in L/mol
					//gm_ph[i*nxyz+j] = max(Cg[i]()[j]*gvol[j]/phreeqcVm[j],1e-6);//Cg in fraction and Vm in L/mol
						}	
					} 
				}
			else 
				{
				for (j=0; j<nxyz;j++) { 
					forAll(Cg,i) {
						gm_ph[i*nxyz+j] = max(Cg[i]()[j]*gvol[j]/phreeqcVm[j],1e-16);
						} //Cg in fraction and Vm in L/mol
					} 
				}
			//for (j=0; j<nxyz;j++) {p_ph[j]=p[j]/atmPa;}
			
			Info<<" phqVm[0] "<<phreeqcVm[0]<<" "<<endl;
			
			//auto start = std::chrono::high_resolution_clock::now();
			//################# RUN PHREEQC   ################
			// set saturations using rchange
			
			Info<<"temp ";
			for (j=0; j<nxyz;j++) {rchange[j] *= sw[j];t_ph[j]=T[j];}//Info<<T[j]<<" ";}//Info<<"rch "<<rchange[j]<<endl;}
			freak.setGvol(gvol); // set gas volume in phreeqc
			freak.setWsat(rchange); // rchange for the calculation doamin, with 0 outside, sw saturation
			freak.setC(c_ph);//transfer c_ph to freak : it does not work to send directly to freak.c
			freak.setGm(gm_ph);//transfer gm_ph to freak
			freak.setTemp(t_ph);
			//freak.setP(p_ph);//transfer pressure to freak
			freak.setTstep(runTime.value()-oldTime); //Info<<" this tme "<< runTime.value()<<" old "<<oldTime<<endl;//the calculation time shall include all time since las phreeqc run
			Info << "running phreeqc dt "<<runTime.value()-oldTime<<endl;
			freak.run();
			freak.getSelOutput();
			Info << "phreeqc done "<<endl;
				
			// write to intermediate file, input (for gases)
			if (ph_gcomp>0) {
				std::ofstream inPhq(cur_dir/"phq_input.txt");
				for (j=0; j<nxyz;j++) { inPhq<<j<<" "<<p[j]<<" "<<rchange[j]<<" "<<gvol[j]<<" "; for (i=0; i<ph_gcomp;i++) {inPhq<< gm_ph[i*nxyz+j] <<" ";} inPhq<<"\n"; }
				}
			
			//auto finish = std::chrono::high_resolution_clock::now();
			//std::chrono::duration<double> dt = finish - start;dure = dure+dt.count();
			// transfer back to C but before keep the previous values for outside domain of calculation
			forAll(Cw,i) {Cw[i]() = Cw[i]().prevIter();}
			if (ph_gcomp>0) {forAll(Cg,i) {Cg[i]() = Cg[i]().prevIter();} }
			forAll(Cw,i) // dissolved
				{
					for (j=0; j<nxyz;j++)
						{
						Cw[i]()[ractive[j]] = freak.c[i*nxyz+j];
						//if (j==imin) {Info<<"ic "<<i<<" imin "<<imin<<" c "<<Cw[i]()[imin]<<endl;}
						}
				} 
			
			// gas, read partial pressures (freak.g in atm) set Cg to fraction (freak.g/Cgtot) and set p to sum of Cg
			if (flowType == 4) // multiphase
				{
				for (j=0; j<nxyz;j++) //should consider ractive
					{
					Cgtot = 0; Gmtot = 0;
					forAll(Cg,i) {freak.g[i*nxyz+j] = max(freak.g[i*nxyz+j],0.); Cgtot += freak.g[i*nxyz+j]; Gmtot += freak.gm[i*nxyz+j];}
					forAll(Cg,i) {Cg[i]()[j] = freak.g[i*nxyz+j]/Cgtot;} // /Cgtot/phreeqcVm;}
					phreeqcVm[j] = gvol[j]/Gmtot;
					p[j] = Cgtot*atmPa;
					//sw[j] = freak.wsat[j];
					}
				for (int i=0; i<ph_gcomp;i++){for (j=0;j<3;j++) {Info <<"g_spc "<< i <<" Cg "<< Cg[i]()[j] <<" sw "<<sw[j]<< endl;}}
				}
			else if (ph_gcomp>0) //unsaturated (not diffrent from above now)
				{
				for (j=0; j<nxyz;j++) //should consider ractive
					{
					Cgtot = 0;
					//phreeqcVm[j] = freak.g[j]*gvol[j]/freak.gm[j];
					forAll(Cg,i) {freak.g[i*nxyz+j] = max(freak.g[i*nxyz+j],0.); Cgtot += freak.g[i*nxyz+j];}
					forAll(Cg,i) {Cg[i]()[j] = freak.g[i*nxyz+j]/Cgtot;} // /Cgtot or Cgtot/phreeqcVm;}
					}
				}
	//Info<<"gvol 20 "<<gvol[20]<<" cg 0 19 "<<Cg[0]()[19]<<" cg 0 20 "<<Cg[0]()[20]<<" cg 0 21 "<<Cg[0]()[21]<<endl;
			// write to phq output file			
			if (ph_gcomp>0) {
				std::ofstream outPhq(cur_dir/"phq_output.txt");
				for (j=0; j<nxyz;j++) { outPhq<<j<<" "<<p[j]<<" "<<freak.gvol[j]<<" "<<phreeqcVm[j]<<" "<<freak.p[j]<<" "<<freak.wsat[j]<<" "; 
					for (i=0; i<ph_gcomp;i++) {outPhq<< freak.g[i*nxyz+j] <<" ";} 
					for (i=0; i<ph_gcomp;i++) {outPhq<< freak.gm[i*nxyz+j] << " ";} 
					outPhq<<"\n"; 
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
					//a1 = a1 *   //gm_ph = Cg[i]()*gvol/phreeqcVm
					if (j<5) {Info<< " gm_ph "<< gm_ph[iw*nxyz+j] << " frk "<< freak.gm[iw*nxyz+j] <<" a1 "<<a1<< endl;}
					sw[j] = max(sw_min[j],sw[j] - a1*.01801/eps[j]); 
					}
				} 
				//nb of moles of H2O(g) transformed in water volume (1 mol 18.01 mL at 25°C)
			for (j=0;j<3;j++) {Info <<" new sw "<<sw[j]<< endl;}
				
			oldTime = runTime.value()*1;
		} //end activate reaction
		
		bool ts;
		ts = runTime.write();//oldTime=runTime.value();
		//write species
		if (ts && flowType==4) {phiGr.write();}
		if (ts && activateReaction==1) {
			phiw.write();phig.write();
			std::ofstream outFile(cur_dir/ runTime.timeName() /"Species");
			outFile.unsetf(std::ios::scientific);outFile.precision(6);
			for (const auto &x : freak.spc) outFile << x << "\n";
			/*//start write species file in columns
			j=0;
			for (i=0;i<cells.size();i++) {
				for (k=0;k<freak.nspc;k++) {
					if (i==ractive[j]) {outFile << freak.spc[i*nxyz+j]<<" ";j+=1;}
					else {outFile << 0.<<" "}
					}
				outFile<<"\n";
			} //end write species file*/
			
			}
		
		Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
			<< "  ClockTime = " << runTime.elapsedClockTime() << " s"
			<< nl << endl;

		Info<< "End\n" << endl;
		double rt = runTime.controlDict().lookupOrDefault("writeInterval",0);
		Info <<" time" << mesh.time().time().value() <<" w intv " <<rt<<endl;
		//if (mesh.time().time().value()>rt/10) {
		//int inpt = cin.get();//}
		istep ++;
	}
    return 0;
}

// ************************************************************************* //
