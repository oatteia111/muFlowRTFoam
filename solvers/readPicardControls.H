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

Global
    readPicardControls

Description
    Read the control parameters used by setDeltaT for solver with Picard's iterations

\*---------------------------------------------------------------------------*/

const dictionary& PicardDict = solutionDict.subDict("Picard");

const scalar tolPicard = PicardDict.lookupOrDefault<scalar>("tolerance", 1.e-3);
//const scalar tolPicard = 1e-3;
const label maxIterPicard = PicardDict.lookupOrDefault<label>("maxIter",10);
//const label maxIterPicard = 12;

//- if dtFromPicard used
const label nIterStability = PicardDict.lookupOrDefault<label>("nIterStability",5);
const label minIterPicard = PicardDict.lookupOrDefault<label>("minIter",1);
//label nIterStability = 6;
//const label minIterPicard = 1;
const scalar dTFactIncrease = runTime.controlDict().lookupOrDefault<scalar>("dTFactIncrease",1.5);
const scalar dTFactDecrease = runTime.controlDict().lookupOrDefault<scalar>("dTFactDecrease",0.5);


//- if dtFromtheta used
const scalar dswmax = runTime.controlDict().lookupOrDefault("dswmax", 5e-3);
scalar dswdTmax = dswmax/runTime.deltaTValue();
scalar dsw = 0;
scalar dsw_avg = 0;

if (adjustTimeStep)
{
    if (timeStepControl == "Picard")
    {
        Info << "Timestepping is based on Picard's iteration with :"
            << nl << "{"
            << nl << "    dTFactIncrease = " << dTFactIncrease
            << nl << "    dTFactDecrease = " << dTFactDecrease
            << nl << "    minIterPicard  = " <<  minIterPicard
            << nl << "    nIterStability = " << nIterStability
            << nl << "}"
            << endl;
    }
    else if (timeStepControl == "dswmax")
    {
        Info << "Timestepping is based on theta variation with :"
            << nl << "{"
            << nl << "    dswmax   = " << dswmax
            << nl << "}"
            << endl;
    }
}
else
{
    Info << "Timestep is fixed";
}

// ************************************************************************* //
