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
#include "fvOptions.H"
#include "fvPatchFieldMapper.H"
#include "simpleControl.H"
// #include "harmonic.H" work only for scalars (wanted for K here)

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

	#include "createFvOptions.H"
	
	dimensionedScalar sT = runTime.startTime();
	dimensionedScalar eT = runTime.endTime();
	dimensionedScalar dT = mesh.time().deltaTValue();
	
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
	//void Time::readDict();
	//runTime.setDeltaT(dT); // 12/3/21 
	//mesh.time().setControls();//readDict();
	//Foam::Time runTime1(Foam::Time::controlDictName, args); // not allowed

	while (runTime.run())
    {
		#include "setDeltaT.H"
		runTime++;//nbt += 1;
		Info << "time = " << runTime.timeName() <<  "  deltaT = " <<  runTime.deltaTValue() << endl;
		// *********** here provide change of density and viscosity if required
		
		//***********************  solve transient flow   *******************************
		if (flowType>0) {
			#include "flow.H"
			}
		bool ts;
		ts = runTime.write();
		//write species
		if (ts) { phiw.write();}
		
		Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
			<< "  ClockTime = " << runTime.elapsedClockTime() << " s"
			<< nl << endl;

		Info<< "End\n" << endl;
		//Info <<"phreeqc time" << dure << endl;
	}
    return 0;
}

// ************************************************************************* //
