Uw.correctBoundaryConditions();
forAll(mesh.boundary(),patchi)
{
	if (isA< fixedValueFvPatchField<vector> >(Uw.boundaryField()[patchi]))
		{ phiw.boundaryFieldRef()[patchi] = Uw.boundaryField()[patchi] & mesh.Sf().boundaryField()[patchi]; }
}

if (flowType == 3) // unsaturated
{
	Info <<"in unsat "<<endl;
   if (timeStepControl == "Picard")
   {
	  iterPicard = 0; resPicard = 1000.;
		while (resPicard > tolPicard)
		{
			iterPicard++;
			#include "hpEqn.H"
			#include "updateThetaFields.H"
		resPicard = gMax((mag(hp-hp.prevIter()))->internalField());
		Info << "Picard iteration : " << iterPicard << " residual = " << resPicard << endl; 
			if (iterPicard == maxIterPicard)
			{
				Warning() <<  " Max iteration reached in Picard loop" << endl;
				break;
			}
		}
    }
	else {
		#include "hpEqn.H"
		#include "updateThetaFields.H"
	volScalarField dsw_tmp = mag(sw-sw.oldTime());
	dsw = gMax(dsw_tmp);
	dswdTmax = dsw/runTime.deltaTValue();Info<<"dswtmp 30 "<<dsw_tmp[30]<<" dsw "<<dsw<<" dswdtmax "<<dswdTmax<<endl;
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
}
else if (flowType != 0) // flowType=0 no flow
{
   Info <<"in simple flow "<<endl;
   #include "hEqn.H"
}
Uw = fvc::reconstruct(phiw); // This is a true darcy velocity (flux divided by surface even in unconfined)
