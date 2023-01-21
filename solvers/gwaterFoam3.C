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
#include <chrono>  // for high_resolution_clock
#include "fvCFD.H"
#include "incompressiblePhase.H"
#include "inletOutletFvPatchField.H"
#include "volFields.H"
#include "fvPatchFieldMapper.H"
#include "fvOptions.H"
#include "simpleControl.H"

#include "phreeqc/initPhreeqc.H"
std::vector<double> a(12,0.);
std::vector<double> c_ph,g_ph,gvol;
int i,j;
my_phq freak; //need to be here to be availabel for every chem condition

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;

int main(int argc, char *argv[])
{
	//init openFoam
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "readGravitationalAcceleration.H"
    #include "createFields.H"
	#include "readPicardControls.H"
    simpleControl simple(mesh);
	#include "createThetaFields.H"
	#include "create2phaseFields.H"
	
	// data input
	std::ifstream inputRactive{cur_dir+"/constant/options/ractive" }; // version 0 shall contain 0 for inactive and 1 for active reaction cell
	std::vector<int>  ractive = {std::istream_iterator<int>{inputRactive}, {}};
	std::ifstream inputInit{cur_dir/"phqinit.txt"};
	std::vector<int> ph_init{std::istream_iterator<int>{inputInit}, {}}; //for (int i=0; i<7;i++){Info << "init nb "<< ph_data[i] << endl;}
	nxyz=ph_init[0];ph_ncomp=ph_init[1];ph_nsolu=ph_init[3]; //Info << "start phreeqc " << endl;
	freak.setDB(cur_dir/"phreeqc.dat");
	//make a first init to calculate the inital solutions
	freak.setData(ph_init);
	std::vector<double> poro(nxyz,0.25);
	for (j=0; j<nxyz;j++){poro[j]=eps[j];}
	freak.setPoro(poro);
	freak.setChemFile(cur_dir/"initChem.pqi"); //Info << "initCh read " << endl;
	freak.init(); //Info << "nxyz " << nxyz << endl;
	
	//writes the inital solutions to file
	std::ofstream outFile(cur_dir/"constant/options/solutions");
	std::vector<float> solu_conc(ph_nsolu*ph_ncomp,0.);Info << "nsolu ncomp "<<ph_nsolu << " "<< ph_ncomp <<endl;
	for (i=0;i<ph_nsolu;i++) // solu number
		{ for (j=0;j<ph_ncomp;j++) // component number
			{
			float a = freak.c[j*ph_nsolu+i];
			outFile << a << "\n"; 
			solu_conc[i*ph_ncomp+j] = a;
			} 
		}
	outFile.close();
	
	//make the true initializaiton for the ractive domain
	std::ifstream inputData{cur_dir/"phqfoam.txt"};
	std::vector<int> ph_data{std::istream_iterator<int>{inputData}, {}}; //for (int i=0; i<7;i++){Info << "init nb "<< ph_data[i] << endl;}
	nxyz=ph_data[0];ph_ncomp=ph_data[1];ph_gcomp=ph_data[2];ph_nsolu=ph_data[3]; //!!! nxyz here is inside the ractive part
	//initiate poro
	poro.resize(nxyz,0.25);
	for (j=0; j<nxyz;j++){poro[j]=eps[j];}
	freak.setPoro(poro);
	gvol.resize(nxyz,1e-8);
	freak.setGvol(gvol); // set gas volume in phreqc
	freak.setData(ph_data); //Info << "db read " << endl;
	freak.init(); //Info << "nxyz " << nxyz << endl;
		
	#include "createFvOptions.H"

	//build the c_ph fields and get conc from phreeqc
	c_ph.resize(nxyz*ph_ncomp,0);
	for (i=0; i<ph_ncomp;i++)
		for (int j=0;j<nxyz;j++) 
			{c_ph[i*nxyz+j] = freak.c[i*nxyz+j];}

	//build the Cwi fields and attribute them concentrations
	#include "transport/createCwiFields.H"
	int icnt = 0;
	forAll(Cw,i) {
		for (j=0; j<nxyz;j++){
			Cw[i]()[ractive[j]] = freak.c[icnt];icnt ++;
			if (j==100) {Info << "sp "<<i << " cell "<< j << " " << Cw[i]()[j] << endl;}
			} // transfer freak.c to Cw
		} 
	//for (i=0; i<ph_ncomp;i++){Info << i << " C[4] "<< Cw[i]()[4] << endl;}
	/*forAll(Cw,i) {
		string s = cur_dir/"0/Cwin";s+= std::to_string(i);
		std::ofstream outFile(s);
		for (const auto &x : Cw[i]()) outFile << x << "\n";
		outFile.flush();outFile.close();
	}*/
	
	//volScalarField b("b",unity);
	//b.dimensions().reset(dimless/dimLength);
	dimensionedScalar sT = runTime.startTime();
	dimensionedScalar eT = runTime.endTime();
	dimensionedScalar dT = mesh.time().deltaTValue();
	//const scalar dt1 = runTime.controlDict().lookupOrDefault("deltaT1", 100);

	if ((flowStartSteady==1.)&&(flowType>0.))
		{
		while (simple.loop(runTime))
			{
			#include "hstdEqn.H"
			}
		}
	runTime.setTime(sT,0); // 12/3/21 time value and index
	runTime.setEndTime(eT); // 12/3/21 set end time that was lost during the simple loop in hstdEqn
	runTime.runTimeModifiable();
	runTime.read();
	float oldTime=0;
	
	int istep = 0;int tcnt = 0;
scalar dC=1e-9;scalar dC1 = 1e-9;scalar dtForC = 1;
	while (runTime.run())
    {
		runTime++;
		Info << "time = " << runTime.timeName() <<  "  deltaT = " <<  runTime.deltaTValue() << endl;
		// *********** here provide change of density and viscosity if required
		
		//***********************  solve transient flow   *******************************
		if (flowType>0) {
			#include "flow.H"
			}
		
		//***************  solve Transport & reaction  *************************
		forAll(Cw,i) {Cw[i]().storePrevIter();}
		#include "transport/CwiEqn.H"
		#include "transport/setDeltaTtrsp.H"
		dC1 = dC*.999;Info<<"dC1  "<<dC1<<endl;
		
		// find where the transported have changed to calculate only there
		icnt = 0;
		std::vector<double> rchange(nxyz,0.);Info <<"ractive size "<<ractive.size()<<endl;
		for (i=4; i<ph_ncomp; i++)
			{
			for (j=0; j<nxyz;j++)
				{
				if (abs(c_ph[i*nxyz+j]-Cw[i]()[ractive[j]])/(c_ph[i*nxyz+j]+1e-20)>1e-5 && sw[ractive[j]]>0.001) {rchange[j] = max(rchange[j],1.);}
				}
			} 
		double ichange = 0;
		for (j=0; j<nxyz;j++) {ichange+=rchange[j];}
		Info <<"nb of changing cells "<<ichange<<endl;
		
		tcnt++;
		if (tcnt>9) {tcnt=0;}
		if (activateReaction==1 && tcnt==0)
		{
			// transfer to phreeqc and run chemistry
			int icnt = 0;
			for (j=0; j<nxyz;j++) 
				{
				//Info << ractive[j] ;
				forAll(Cw,i) 
					{
					c_ph[i*nxyz+j] = Cw[i]()[ractive[j]];
					} 
				//Info << " "<< endl;
				} // transfer C to c_ph
			
			//auto start = std::chrono::high_resolution_clock::now();
			
				freak.setC(c_ph);//transfer c_ph to freak : it does not work to send directly to freak.c
				freak.setWsat(rchange);
				freak.setTstep(runTime.deltaT().value()-oldTime);
				Info << "running phreeqc "<<endl;
				freak.run();
				freak.getSelOutput();
				Info << "phreeqc done "<<endl;
			
			//auto finish = std::chrono::high_resolution_clock::now();
			//std::chrono::duration<double> dt = finish - start;dure = dure+dt.count();
			
			// transfer back to C but before reset all to 0 (for outside domain of calculation
			forAll(Cw,i) {Cw[i]() = Cw[i]().prevIter();}
			forAll(Cw,i) {for (j=0; j<nxyz;j++)
				{
				Cw[i]()[ractive[j]] = freak.c[i*nxyz+j];
				//if (ractive[j]%ncell_lay==6032) {Info << "i "<<i<<" j "<<j<<" cell "<<ractive[j]<<" c_ph "<< c_ph[i*nxyz+j] <<" freakc "<< freak.c[i*nxyz+j] <<endl;}
				}} 
			//for (int i=0; i<ph_ncomp;i++){Info << i <<"C 118, 1732 "<< Cw[i]()[118] << "  " << Cw[i]()[1732] << endl;}
			oldTime = runTime.value();
		}
		
		bool ts;
		ts = runTime.write();
		//write species
		if (ts ) {
			phiw.write();
			std::ofstream outFile(cur_dir/ runTime.timeName() /"Species");
			for (const auto &x : freak.spc) outFile << x << "\n";
		}
		
		Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
			<< "  ClockTime = " << runTime.elapsedClockTime() << " s"
			<< nl << endl;

		Info<< "End\n" << endl;
		//Info <<"phreeqc time" << dure << endl;
		istep ++;
	}
    return 0;
}

// ************************************************************************* //
