 //- Choose the method to compute time step
/*const bool adjustTimeStep =
    runTime.controlDict().lookupOrDefault("adjustTimeStep", false);
const word timeStepControl =
    runTime.controlDict().lookupOrDefault<word>("timeStepControl", "Picard");*/

//- Time control with C variation
/*scalar dCmax_abs =
    runTime.controlDict().lookupOrDefault<scalar>("dCmax_abs",1e-6);
scalar dCmax_rel =
    runTime.controlDict().lookupOrDefault<scalar>("dCmax_rel",1e-6);*/
//Info<< "dCm_a " << dCmax_abs << " dCm_r " << dCmax_rel << endl;

//- C-variation control
//dimensionedScalar dCmax(transportProperties.lookupOrDefault("dCmax",dimensionedScalar("dCmax",dimMass/dimVol,0.01)));
//const scalar dCresidual = runTime.controlDict().lookupOrDefault("dCresidual", 1e-3);
scalar maxDeltaT = runTime.controlDict().lookupOrDefault("maxDeltaT", GREAT);
const scalar dCmax = runTime.controlDict().lookupOrDefault("dCmax", 0.01);
Info<< "dCmax " <<dCmax<<" dC1 "<< dC <<" deltaTValue " << runTime.deltaTValue()<<endl;

//scalar dCdTmax = max(dCmax_rel*gMax(Cw),dCmax_abs)/runTime.deltaTValue();
//dtForC =  max(dCmax*gMax(Cw),dCresidual)/(dCdTmax+SMALL);
//if (dC1>dCmax) {dtForC = runTime.deltaTValue();}
//else {dtForC *=1.2;}
dtForC = dCmax/(max(dC1,0)+SMALL)*runTime.deltaTValue();
Info<< "dtForC " << dtForC << endl; 

scalar newDeltaT = min(dtForC, 1.2*runTime.deltaTValue());
//Info<< "newDeltaT " << newDeltaT << " maxDeltaT " << maxDeltaT << endl;

runTime.setDeltaT (min (newDeltaT,maxDeltaT) );

//Info<< "newDeltaT " << newDeltaT << " maxDeltaT " << maxDeltaT << endl;

// ************************************************************************* //