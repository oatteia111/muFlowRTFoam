scalar CFLUse = -1;
scalar maxDeltaTFact = -1;
scalar CoNum = 0.0;
scalar meanCoNum = 0.0;

const scalar maxCo = runTime.controlDict().lookupOrDefault("maxCo", 0.75);
const scalar maxDeltaT = runTime.controlDict().lookupOrDefault("maxDeltaT", 1000);

if (mesh.nInternalFaces())
{
    scalarField sumPhiw
        (fvc::surfaceSum(mag(phiw))().internalField());

    CoNum = 0.5*gMax(sumPhiw/mesh.V().field())*runTime.deltaTValue();

    meanCoNum = 
        0.5*(gSum(sumPhiw)/gSum(mesh.V().field()))*runTime.deltaTValue();
}

Info<< "Courant Number mean: " << meanCoNum << " max: " << CoNum << endl;

CFLUse = CoNum;
maxDeltaTFact = maxCo/(CFLUse + SMALL);

scalar deltaTFact = min(min(maxDeltaTFact, 1.0 + 0.1*maxDeltaTFact), 1.2);

runTime.setDeltaT
(
    min(deltaTFact*runTime.deltaTValue(),maxDeltaT)
);
//Info<< "deltaT " << runTime.deltaTValue() << endl;
//Info<< "end setDeltaT" << endl;

// ************************************************************************* //
