scalar CFLUse = -1;
scalar deltaTFact = 1;
scalar maxDeltaTFact = -1;
scalar CoNum = 0.0;
scalar meanCoNum = 0.0;
scalar deltaT = 0.;
const scalar maxCo = runTime.controlDict().lookupOrDefault("maxCo", 0.75);
scalar maxDeltaT = runTime.controlDict().lookupOrDefault("maxDeltaT", GREAT);
scalar dswmax = runTime.controlDict().lookupOrDefault("dswmax", 0.02);

Uw.correctBoundaryConditions();
forAll(mesh.boundary(),patchi)
{
	if (isA< fixedValueFvPatchField<vector> >(Uw.boundaryField()[patchi]))
		{ phiw.boundaryFieldRef()[patchi] = Uw.boundaryField()[patchi] & mesh.Sf().boundaryField()[patchi]; }
}

if (flowType == 3) // unsaturated
{
   if (timeStepControl == "Picard")
   {
	   Info<<"in picard "<<endl;
	iterPicard = 0; resPicard = 1000.;
	while (resPicard > tolPicard)
		{
		iterPicard++;
		#include "hpEqn.H"
		#include "updateThetaFields.H"
		//resPicard = gMax((mag(hp-hp.prevIter()))->internalField());
		Info << "Picard iteration : " << iterPicard << " residual = " << resPicard << endl; 
		if (iterPicard == maxIterPicard)
			{ Warning() <<  " Max iteration reached in Picard loop" << endl; 
			break; }
		}
	if (iterPicard == maxIterPicard)
		{deltaTFact = 0.8;}
	else if (iterPicard <= nIterStability)
		{deltaTFact = 1.2;}
	runTime.setDeltaT (min(deltaTFact*runTime.deltaTValue(),maxDeltaT));

    }
	
	else 
	{
		#include "hpEqn.H"
		#include "updateThetaFields.H"
		//Info<<"update done "<<endl;
		volScalarField::Internal dsw_tmp = mag(sw-sw.oldTime());
		dsw= gMax(dsw_tmp)/average(sw).value(); // cannot do directly linne aboe and that one
		//dswdTmax = dsw/runTime.deltaTValue();
		scalar ratio = dswmax/(dsw+SMALL); Info<<"dsw "<<dsw<<" dswmax "<<dswmax<<" ratio "<<ratio<<endl;
		/* //not realy able to rewind time
		if (ratio<0.5 && runTime.timeIndex()>2) {
			hp = hp.oldTime();
			Info<<"time "<<runTime.timeOutputValue()<<" dt "<<runTime.deltaTValue()<<" index "<<runTime.timeIndex()<<endl;
			runTime.setTime(runTime.timeOutputValue()-runTime.deltaTValue(),runTime.timeIndex());
			//deltaTFact = min(ratio,50.);//runTime.setDeltaT(runTime.deltaT()/min(ratio,50.));
			runTime.setDeltaT(min(max(min(ratio,50.)*runTime.deltaTValue(),1.2),maxDeltaT));
			}
		else {
		*/
			deltaTFact = max(min(ratio,1.2),0.1);
			runTime.setDeltaT(min(deltaTFact*runTime.deltaTValue(),maxDeltaT));
			//runTime.setTime(runTime.timeOutputValue()+runTime.deltaTValue(),runTime.timeIndex());
			//}
	//runTime.setDeltaT(min(max(deltaTFact*runTime.deltaTValue(),1),maxDeltaT));
    }
}

else if (flowType == 4) // 2 phase
{
   Info <<"in 2phase "<<endl;
	Ug.correctBoundaryConditions();
	forAll(mesh.boundary(),patchi)
	{
		if (isA< fixedValueFvPatchField<vector> >(Ug.boundaryField()[patchi]))
			{ phig.boundaryFieldRef()[patchi] = Ug.boundaryField()[patchi] & mesh.Sf().boundaryField()[patchi]; }
	} 
    #include "pSEqn.H"
	#include "CoatsNo.H"
	deltaTFact = min(min(maxDeltaTFact, 1.0 + 0.1*maxDeltaTFact), 1.1);
	runTime.setDeltaT (min(deltaTFact*runTime.deltaTValue(),maxDeltaT));
}
else if (flowType != 0) // flowType=0 no flow
{
   Info <<"in simple flow "<<endl;
   #include "hEqn.H"
   deltaTFact = 1.2;
   runTime.setDeltaT (min(deltaTFact*runTime.deltaTValue(),maxDeltaT));
}
//Info<<" deltaTFact " <<deltaTFact<<" maxDeltaT "<<maxDeltaT<<" tnext(s) "<<tnext*day<<endl;
//maxDeltaT = min(maxDeltaT,(tnext-runTime.value())*day);

Uw = fvc::reconstruct(phiw); // This is a true darcy velocity (flux divided by surface even in unconfined)
