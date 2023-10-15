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
#include <chrono>

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

#include "nnBase.h"

#include <chrono>  // for high_resolution_clock
#include "fvCFD.H"
#include "incompressiblePhase.H"
#include "inletOutletFvPatchField.H"
#include "volFields.H"
#include "fvPatchFieldMapper.H"
// #include "fvOptions.H"
#include "simpleControl.H"
#include "cellSet.H"
#include "phreeqc/initPhreeqc.H"

std::vector<double> a(12,0.);
std::vector<double> c_ph,gm_ph,p_ph,S_ph,gvol,ractive,solu_conc,gas_conc;
float atmPa=101325.;float vmw,Cgtot,Gmtot;
int i,j,iw;double x;
my_phq freak; //need to be here to be availabel for every chem condition
#include "myFunc.H"
dicFunc fDe_T;
//create NN
my_NN Cwgnn;	   		  

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

std::vector<int> indexC(labelList &cells, std::vector<float> &data)
{
    std::vector<int> c1(cells.size(),0);//Info<<"in index "<<cells.size()<<endl;
	for (i=0; i<cells.size();i++)  // reads the first ncells lines
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

outData getCbuffer(string fname, int itime, int ncell) {
    std::vector<float> data(ncell*4);
	std::ifstream inputData{cur_dir+"/constant/options/"+fname, std::ios::binary}; //
	inputData.seekg(ncell*itime*4*sizeof(float)); //each line is composed of 4 numbers and there are two values at the beginning

    inputData.read(reinterpret_cast<char*>(&data[0]), ncell*4*sizeof(float));
	float time;	
	inputData.read(reinterpret_cast<char*>(&time), sizeof(float));
	outData output;
	output.t = time;std::cout<<"readbin "<<fname <<" "<<itime<<" "<<ncell;
	output.d = data;
    return output;
}

//----------------------------- MAIN ------------------------------

using namespace Foam;
int main(int argc, char *argv[])
{
	namespace chr = std::chrono;
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
	//for NN	
	std::random_device rd;
	std::mt19937 shuf(rd()); // will be used later
	std::vector<double> Cmin,Cmax,cdff,Ymin,Ymax; // Cmin,max for conc, Y for selected species
	std::vector<int> llog; //to store the variable that are log transformed (for nn)
	int nncols,nsel;
	char sep ='/'; //don't understand why here sep is needed as char while in muFlow just work directly
	std::ofstream outTime(cur_dir+sep+"time.txt");

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

	//#################  NN init ################### (must be before full domain as it also initialize cells in phreeqc)
	//Cwgnn.setRunParms({nn_epoc,nn_batch,nn_lr,0.75});Cwgnn.init(); //+1 for time step
	std::ofstream outNNrmse(cur_dir/"NNrmse.txt");
	std::vector<float> nndata; // TORCH only accepts floats ???? (to vbe validated)
	std::vector<float> nntarget;
	std::ofstream outNNdata(cur_dir/"NNdata.txt");
	std::ofstream outNNtarget(cur_dir/"NNtarget.txt");
	if (activateNNchemistry==1) 
		{
		#include "phreeqc/nnFreakFirst.H"
		}
	
	//##################### make the initialization for the full domain in phreeqc : data,poro, gvol
	std::ifstream inputData{cur_dir/"phqfoam.txt"};
	std::vector<int> ph_data{std::istream_iterator<int>{inputData}, {}}; //for (int i=0; i<7;i++){Info << "init nb "<< ph_data[i] << endl;}
	nxyz=ph_data[0];ph_ncomp=ph_data[1];ph_gcomp=ph_data[2];ph_nsolu=ph_data[3]; //!!! nxyz here is inside the ractive part
	freak.setData(ph_data); Info << "nxyz " << nxyz << endl;
	freak.setDB(cur_dir+sep+"phreeqc.dat");
	freak.setChemFile(cur_dir/"initChem.pqi"); //Info << "initCh read " << endl;
	//initiate poro and gas volume
	poro.resize(nxyz,0);
	for (i=0;i<nxyz;i++) {poro[i]=eps[i];}
	std::vector<double> wsat(ph_nsolu,0.9999);
	freak.setPoro(poro);
	freak.setWsat(wsat); // rchange for the calculation doamin, with 0 outside, sw saturation
	//for (i=0;i<nxyz;i++) {p_ph[i]=p[i]/atmPa;}	
	//freak.setP(p_ph);
	freak.init();
	gvol.resize(nxyz,0.01);
	//freak.run();

	
	//##############" build the c_ph and gm_ph fields and get conc from phreeqc (c_ph=freak.c but needed two variables for format questions)
	c_ph.resize(nxyz*ph_ncomp,0);
	for (i=0; i<ph_ncomp;i++)
		for (j=0;j<nxyz;j++) 
			{c_ph[i*nxyz+j] = freak.c[i*nxyz+j];}
	// gases are in atm (freak.g) we start with ideal gas
	gm_ph.resize(nxyz*ph_gcomp,0); // moles of gas
	Info<<"gcomp "<<ph_gcomp<<" g size "<<freak.g.size()<<" frk.c size "<<freak.c.size()<<endl;
	for (i=0; i<ph_gcomp;i++)
		for (j=0;j<nxyz;j++) 
			{
			phreeqcVm[j] = 24.5/(p[j]/atmPa);
			gm_ph[i*nxyz+j] = freak.gm[i*nxyz+j];
			}
	iw = freak.iGwater;
	// also store results of selectedoutput (for nn) and recalc Cmin Cmax (only to extend them if needed)
	if (activateNNchemistry==1) 
		{
		S_ph.resize(nxyz*nsel);freak.getSelOutput();nsel=freak.nselect-1; // S_ph will store all data on surface
		for (j=0;j<nxyz;j++)
			{
			S_ph[j]=freak.spc[j];//pH, we remove pe
			for (i=1;i<nsel;i++) {S_ph[i*nxyz+j]=freak.spc[(i+1)*nxyz+j];}
			}
		for (j=0;j<nxyz;j++)
			{
			for (i=0;i<nncols;i++)
				{
				if (i<ph_ncomp-4) {x = c_ph[i*nxyz+j]; }
				else {x=S_ph[i*nxyz+j]; }
				x = std::max(x,1e-12);
				if (llog[i]==1) {Cmin[i]=std::min(Cmin[i],std::log10(x));Cmax[i]=std::max(Cmax[i],std::log10(x));}
				else {Cmin[i]=std::min(Cmin[i],x);Cmax[i]=std::max(Cmax[i],x);}
				}
			}

		std::cout<<"Cmin "; for (i=0;i<nncols;i++) {std::cout<<Cmin[i]<<" ";} ; std::cout<<"\n";
		std::cout<<"Cmax "; for (i=0;i<nncols;i++) {std::cout<<Cmax[i]<<" ";} ; std::cout<<"\n";
		}
	
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
	int tcnt=0;int flag1step=0;	int istep = 0;
	std::vector<double> rchange(nxyz,0.);
	
	while (runTime.run())
    {
		runTime.read();
		// #include "transport/setDeltaTtrsp.H"
		runTime++;
		Info << "time = " << runTime.timeName() <<  "  deltaT = " <<  runTime.deltaTValue() << endl;
		outTime << mesh.time().value()/86400.<<" ";
		// *********** here provide change of density and viscosity if required
		
		//***********************  solve transient flow   *******************************
		//for (j=0; j<nxyz;j++) {if (j<6) {Info<<"p before flow "<<p[j]/atmPa<<" sw "<<sw[j]<<endl;}}
		if (flowType>0) {
			auto t1 = chr::high_resolution_clock::now();
			#include "flow.H"
			auto t2 = chr::high_resolution_clock::now();
			auto ms_int = chr::duration_cast<chr::milliseconds>(t2 - t1);
			outTime << ms_int.count()<< " ";
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
				auto t1 = chr::high_resolution_clock::now();
				forAll(Cw,i) {Cw[i]().storePrevIter();} // for cells outside calculation
				#include "transport/CwiEqn.H"
				auto t2 = chr::high_resolution_clock::now();
				auto ms_int = chr::duration_cast<chr::milliseconds>(t2 - t1);
				outTime << ms_int.count()<< " ";
				
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
			auto t1 = chr::high_resolution_clock::now();
			#include "chem_nn.h"
			auto t2 = chr::high_resolution_clock::now();
			auto ms_int = chr::duration_cast<chr::milliseconds>(t2 - t1);
			outTime << ms_int.count()<< "\n";std::flush(outTime);
			oldTime = runTime.value()*1;
		} //end activate reaction
		
		bool ts;std::cout<<"flag1 "<<flag1step<<"\n";
		ts = runTime.write();//oldTime=runTime.value();
		//write species
		if (ts && flowType==4) {phiGr.write();}
		if (ts && activateReaction==1) {
			phiw.write();phig.write();flag1step=1;
			std::ofstream outFile(cur_dir/ runTime.timeName() /"Species");
			outFile.unsetf(std::ios::scientific);outFile.precision(6);
			for (const auto &x : freak.spc) outFile << x << "\n";
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
