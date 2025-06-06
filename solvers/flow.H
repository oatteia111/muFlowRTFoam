scalar CFLUse = -1;
scalar deltaTFact = 1;
scalar maxDeltaTFact = -1;
//scalar CoNum = 0.0;
//scalar meanCoNum = 0.0;
//const scalar maxCo = runTime.controlDict().lookupOrDefault("maxCo", 0.75);
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
	// this loop is imitated from Orgogozo Richards
	deltaTFact = 1.;
	hp1 = hp;
	while (resPicard > tolPicard)
		{
		iterPicard++;
		hp = hp_tmp; 
		#include "flow/updateThetaFields.H"
		#include "flow/hpEqn.H"
		resPicard = gMax((mag(hp-hp_tmp))->internalField());
		hp_tmp = hp;
		if (iterPicard == maxIterPicard)
			{ runTime.setTime(oldTime,tstep); break; }
		}  // end picard iter
	Info << "Picard nb iterations : "<<iterPicard<<endl;
	if (iterPicard == maxIterPicard) {deltaTFact = 0.2;hp_tmp=hp1;hp=hp1;rewind=1;}
	else if (iterPicard <= nIterStability) {deltaTFact = 1.2;}
	newDeltaT = min(deltaTFact*runTime.deltaTValue(),maxDeltaT);
    }
	
	else 
	{
		#include "flow/hpEqn.H"
		#include "flow/updateThetaFields.H"
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
			newDeltaT = min(deltaTFact*runTime.deltaTValue(),maxDeltaT);
			//runTime.setTime(runTime.timeOutputValue()+runTime.deltaTValue(),runTime.timeIndex());
			//}
	//runTime.setDeltaT(min(max(deltaTFact*runTime.deltaTValue(),1),maxDeltaT));
    }
}

else if (flowType == 4) // 2 phase
{
   Info <<"in 2phase "<<endl;
	/*Ug.correctBoundaryConditions();
	forAll(mesh.boundary(),patchi)
	{
		if (isA< fixedValueFvPatchField<vector> >(Ug.boundaryField()[patchi]))
			{ phig.boundaryFieldRef()[patchi] = Ug.boundaryField()[patchi] & mesh.Sf().boundaryField()[patchi]; }
	} */
	//Info<<p<<endl;
    if (timeStepControl == "Picard")
    {
	Info<<"in picard "<<endl;
	iterPicard = 0; resPicard = 1000.;
	// this loop is imitated from Orgogozo Richards
	deltaTFact = 1.;
	p1 = p;
	while (resPicard > tolPicard)
		{
		iterPicard++;
		p = p_tmp; 
		#include "flow/pSEqn.H"
		Info<<" end pS ,p_tmp(0) "<<p_tmp[0]<<endl;
		resPicard = gMax((mag(p-p_tmp)/atmPa*1e3)->internalField());Info<<"res Pic "<<resPicard<<endl;
		p_tmp = p;
		if (iterPicard == maxIterPicard){ break; }
		}  // end picard iter
	Info << "Picard nb iterations : "<<iterPicard<<endl;
	if (iterPicard == maxIterPicard) {deltaTFact = 0.2;p_tmp=p1;p=p1;rewind=1;}
	else if (iterPicard <= nIterStability) {deltaTFact = 1.2;}
	newDeltaT = min(deltaTFact*runTime.deltaTValue(),maxDeltaT);
    }
	
	else 
	{
    #include "flow/pSEqn.H"
	//#include "CoatsNo.H"
	volScalarField::Internal dsw_tmp = mag(sw-sw.oldTime());
	for (j=0; j<nxyz;j++) {if (bcPS[j]==0) {dsw_tmp[j]=0;} }
	dsw= gMax(dsw_tmp)/average(sw).value(); // cannot do directly linne aboe and that one
	scalar ratio = dswmax/(dsw+SMALL);
	deltaTFact = max(min(ratio,1.2),0.1);
	newDeltaT = min(deltaTFact*runTime.deltaTValue(),newDeltaT);
	}
	Info<<" dswmax "<<dswmax<<" ndetlaT "<<newDeltaT<<endl;
}
else if (flowType != 0) // flowType=0 no flow
{
   Info <<"in simple flow "<<endl;
   #include "flow/hEqn.H"
	newDeltaT = min(maxDeltaT, 1.5*runTime.deltaTValue());
}
Info<<" deltaTFact " <<deltaTFact<<" maxDeltaT "<<maxDeltaT<<" newdeltat "<<newDeltaT<<endl;

Uw = fvc::reconstruct(phiw); // This is a true darcy velocity (flux divided by surface even in unconfined)
