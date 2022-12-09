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
// #include "sourceEventFile.H"
// #include "outputEventFile.H"
// #include "patchEventFile.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;

int main(int argc, char *argv[])
{
	//init openFoam
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "readGravitationalAcceleration.H"
    #include "createFieldsT.H"
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	while (runTime.run())
    {
		#include "setDeltaT.H"
		runTime++;
		Info << "time = " << runTime.timeName() << nl << endl;

		Info<< "deltaT = " <<  runTime.deltaTValue() << endl;
		//here provide change of density and viscosity if required
		
		//***********************  solve flow   *******************************
    fvScalarMatrix hEqn
        (
            -fvm::laplacian(Mf,h) //+ fvc::div(phiG)
        );

		hEqn.solve();

		phiw = hEqn.flux() ;//+ phiG;

		Uw = fvc::reconstruct(phiw);
		Uw.correctBoundaryConditions();

		//***************  solve Transport & reaction if needed *************************
		if (runTime.write()) {
			h.write();phiw.write();Uw.write();
		}
		Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
			<< "  ClockTime = " << runTime.elapsedClockTime() << " s"
			<< nl << endl;

		Info<< "End\n" << endl;
		//Info <<"phreeqc time" << dure << endl;
	}
    return 0;
}

// ************************************************************************* //
