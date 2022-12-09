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
#include <chrono>  // for high_resolution_clock
#include "fvCFD.H"
#include "incompressiblePhase.H"
#include "inletOutletFvPatchField.H"
#include "volFields.H"
#include "fvPatchFieldMapper.H"
#include "fvOptions.H"
#include "simpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//////////////////// find local dir
#include <stdlib.h>
#include <unistd.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <iterator>
#include <sstream>
#include <string>
#define GetCurrentDir getcwd

std::string get_current_dir() {
   char buff[FILENAME_MAX]; //create string buffer to hold path
   GetCurrentDir( buff, FILENAME_MAX );
   std::string current_working_dir(buff);
   return current_working_dir;
}
string cur_dir = get_current_dir();

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

	//build the C0 fields
	std::ifstream inputCactive{cur_dir+"/constant/options/cactive" }; // version 0 shall contain 0 for inactive and 1 for active reaction cell
	std::vector<int>  cactive = {std::istream_iterator<int>{inputCactive}, {}};
	volScalarField cmult("cmult",unity*0.);
	if (cactive.size()>0) 
		{for (j=0; j<cactive.size();j++) {cmult[cactive[j]]=1;} }
	//Info<< " Cmult "<<cmult[203870]<<endl;
	#include "transport/createCwFields.H"

	#include "createFvOptions.H"
	
	dimensionedScalar sT = runTime.startTime();
	dimensionedScalar eT = runTime.endTime();
	dimensionedScalar dT = mesh.time().deltaTValue();

	if (flowStartSteady==1 && flowType>0)
		{
		#include "hstdEqn.H"
		Info <<"final t H "<<runTime.endTime().value()<<endl;
		std::ofstream outFile(cur_dir/"endSteadyF");
		outFile << runTime.endTime().value() << "\n"; outFile.close();
		phiw.write();h.write();
		}

	if (transportSteady==1)
		{
		#include "transport/CwEqnStd.H"
		Info <<"final t T "<<runTime.endTime().value()<<endl;
		std::ofstream outFile(cur_dir/"endSteadyT");
		outFile << runTime.endTime().value() << "\n"; outFile.close();
		}
	runTime.setTime(sT,0); // 12/3/21 time value and index
	runTime.setEndTime(eT); // 12/3/21 set end time that was lost during the simple loop in hstdEqn
	runTime.runTimeModifiable();
	runTime.read();

//dimensionedScalar dC=Cw[0];//("dC",dimMass/dimVol,1e-9);
//dimensionedScalar dC1=Cw[0]; //("dC1",dimMass/dimVol,1e-9);Info<<"dC1 "<<dC1<<endl;
scalar dC=1e-9;scalar dC1 = 1e-9;scalar dtForC = 1;
	
	if (transportSteady==0)
	{
		while (runTime.run())
		{
			Info<<"dC1 before dt "<<dC1<<endl;
			#include "transport/setDeltaTtrsp.H"
			runTime.read();
			runTime++;
			Info << "time = " << runTime.timeName() <<  "  deltaT = " <<  runTime.deltaTValue() << endl;
			// *********** here provide change of density and viscosity if required
			
			//***********************  solve transient flow   *******************************
			if (flowType>0.) {
				#include "flow.H"
				}
			//***************  solve Transport  *************************
			#include "transport/CwEqn.H"
			//Info<< " Cw infilt "<<Cw[203870]<<endl;
			if (cactive.size()>0) {Cw = Cw*cmult;} //reset to 0 outside the calculation domain
			//Info<< " Cw infilt "<<Cw[203870]<<endl;
			bool ts;
			ts = runTime.write();
			//write species
			if (ts) {
				Cw.write();
				if (flowType>0) {phiw.write();h.write();}
				else {h.read();phiw.read();}
			}
			
			dC1 = dC*.999;Info<<"dC1  "<<dC1<<endl;
			Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
				<< "  ClockTime = " << runTime.elapsedClockTime() << " s"
				<< nl << endl;

			Info<< "End\n" << endl;
			//Info <<"phreeqc time" << dure << endl;
		}
	}
    return 0;
}

// ************************************************************************* //
