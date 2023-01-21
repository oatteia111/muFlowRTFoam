scalar CFLUse = -1;
scalar deltaTFact = 1;
scalar maxDeltaTFact = -1;
scalar CoNum = 0.0;
scalar meanCoNum = 0.0;

const scalar maxCo = runTime.controlDict().lookupOrDefault("maxCo", 0.75);
scalar maxDeltaT = runTime.controlDict().lookupOrDefault("maxDeltaT", GREAT);
Info<<"write tstep "<<runTime.controlDict().lookupOrDefault("timeStep",0)<<endl;

if (timeStepControl == "Picard")
    {
        scalar deltaTFact=1;
        if(iterPicard == maxIterPicard)
        {
            deltaTFact = dTFactDecrease;
            iterStability=0;
        }
        else if (iterPicard <= minIterPicard)
        {
            iterStability++;
        }

        if (iterStability == nIterStability)
        {
            deltaTFact = dTFactIncrease;
            iterStability = 0;
        }
    }
else if(timeStepControl == "dswmax")
    {
        scalar dtForSw = dswmax / (dswdTmax+SMALL);

        runTime.setDeltaT
            (
                min(
                    dtForSw,
                    min(
                        1.2*runTime.deltaTValue(),
                        maxDeltaT
                    )
                )
            );
    }

else if (flowType==4 && CFLint==0)
	{
		#include "CoatsNo.H"
		deltaTFact = min(min(maxDeltaTFact, 1.0 + 0.1*maxDeltaTFact), 1.05);
		Info<< "deltaTfact " << deltaTFact << endl;
	}
else //classical courant No
    {
	if (mesh.nInternalFaces())
	{
		scalarField sumPhiw
			(fvc::surfaceSum(mag(phiw))().internalField());

		CoNum = 0.5*gMax(sumPhiw/mesh.V().field())*runTime.deltaTValue();

		meanCoNum = 
			0.5*(gSum(sumPhiw)/gSum(mesh.V().field()))*runTime.deltaTValue();
	}

	CFLUse = CoNum;
	maxDeltaTFact = maxCo/(CFLUse + SMALL);
	Info<< "Courant Number mean: " << meanCoNum << " max: " << CoNum <<" maxDeltaTFact "<<maxDeltaTFact<< endl;

	deltaTFact = min(min(maxDeltaTFact, 1.0 + 0.1*maxDeltaTFact), 1.2);
	
	}

		
runTime.setDeltaT
(
	min(deltaTFact*runTime.deltaTValue(),maxDeltaT)
);

//Info<< "end setDeltaT" << endl;

// ************************************************************************* //
