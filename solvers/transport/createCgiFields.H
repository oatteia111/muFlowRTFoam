dimensionedScalar Dg0(transportProperties.lookupOrDefault("Dg0",dimensionedScalar("Dg0",dimArea/dimTime,1e-6)));
scalar diffusionEqn(transportProperties.lookupOrDefault<scalar>("diffusionEqn",0));

Info<< "Reading field Deffg if present\n" << endl;
volTensorField Deffg
(
    IOobject
    ("Deffg",runTime.timeName(),mesh,IOobject::READ_IF_PRESENT,IOobject::NO_WRITE),
    mesh,
	dimensionedTensor("",dimArea/dimTime,Zero) 
);

Info <<"in Cgfield nxyz "<<nxyz<< " g comp "<<ph_gcomp << endl;

// Initialize BCs list for Cc to zero-gradient
/*
wordList bcTypes
	(
	Kh.boundaryField().size(),
	zeroGradientFvPatchScalarField::typeName // works
	);
	
List<word> bcNames;bcNames.append("bc3");
for (int i=0;i<bcNames.size();i++)
	{
	const label inletID = mesh.boundaryMesh().findPatchID(bcNames[i]);
	bcTypes[inletID] = fixedValueFvPatchScalarField::typeName;
	}
*/
List <autoPtr <volScalarField> > Cg(ph_gcomp);
forAll (Cg, i)
{
		word nameCi ("Cg" + name(i));
		Cg[i].reset
		(
			new volScalarField
			(
			 IOobject(nameCi,runTime.timeName(),mesh,IOobject::READ_IF_PRESENT,IOobject::AUTO_WRITE),
			 mesh,
			 dimensionedScalar("",dimMass/dimVol,0.), //freak.solu_conc[i]
			 bcTypes
			)
		);
}

Info <<"end Cgifields " << endl;