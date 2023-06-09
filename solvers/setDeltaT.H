scalar CFLUse = -1;
scalar deltaTFact = 1;
scalar maxDeltaTFact = -1;
scalar CoNum = 0.0;
scalar meanCoNum = 0.0;
scalar deltaT = 0.;
const scalar maxCo = runTime.controlDict().lookupOrDefault("maxCo", 0.75);
scalar maxDeltaT = runTime.controlDict().lookupOrDefault("maxDeltaT", GREAT);
scalar dCmax = runTime.controlDict().lookupOrDefault("dCmax", 0.01);
scalar dTmax = runTime.controlDict().lookupOrDefault("dTmax", 0.01);
dtForC = maxDeltaT;

if (flowType==4 && CFLint==0)
	{
	#include "CoatsNo.H"
	deltaTFact = min(min(maxDeltaTFact, 1.0 + 0.1*maxDeltaTFact), 1.1);
	}
else //classical courant No
    {
	if (mesh.nInternalFaces())
	{
		scalarField sumPhiw
			(fvc::surfaceSum(mag(phiw))().internalField());
		CoNum = 0.5*gMax(sumPhiw/mesh.V().field())*runTime.deltaTValue();
		meanCoNum = 0.5*(gSum(sumPhiw)/gSum(mesh.V().field()))*runTime.deltaTValue();
	}
	CFLUse = CoNum;
	maxDeltaTFact = maxCo/(CFLUse + SMALL);
	deltaTFact = min(min(maxDeltaTFact, 1.0 + 0.1*maxDeltaTFact), 1.2);
	}

deltaT = min(deltaTFact*runTime.deltaTValue(),maxDeltaT);


if (activateTransport==1) 
	{dtForC = dCmax/(max(dC,0)+SMALL)*runTime.deltaTValue();} // we use dtForC for both temperature and C

if (activateThermal==1) 
	{dtForC = min(dtForC,dTmax/(max(dT,0)+SMALL)*runTime.deltaTValue());}

deltaT = min(deltaT, dtForC);
runTime.setDeltaT (deltaT);
// ************************************************************************* //
