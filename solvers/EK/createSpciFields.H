dimensionedScalar alphaL(transportProperties.lookupOrDefault("alphaL",dimensionedScalar("alphaL",dimLength,0.0)));
dimensionedScalar alphaT(transportProperties.lookupOrDefault("alphaT",dimensionedScalar("alphaT",dimLength,0.0)));
dimensionedScalar Dw0(transportProperties.lookupOrDefault("Dw0",dimensionedScalar("Dw0",dimArea/dimTime,1e-10)));

Info<< "Reading field Deffw if present\n" << endl;
volTensorField Deffw
(
    IOobject
    ("Deffw",runTime.timeName(),mesh,IOobject::READ_IF_PRESENT,IOobject::NO_WRITE),
    mesh,
	dimensionedTensor("",dimArea/dimTime,Zero) 
);

// below from https://www.cfd-online.com/Forums/openfoam-programming-development/113261-multiple-volscalarfields-using-list-volscalarfield.html
Info <<"in Cifield nxyz "<<nxyz<< " n comp "<<ph_ncomp << endl;

// Initialize BCs list for Cc to zero-gradient
wordList bcTypes
	(
	Kh.boundaryField().size(),
	zeroGradientFvPatchScalarField::typeName // works
	);
	

List <autoPtr <volScalarField> > ECw(ph_nspc); // here Cw are the water conc of species
forAll (ECw, i)
{
		word nameCi ("ECw" + name(i));
		ECw[i].reset
		(
			new volScalarField
			(
			 IOobject(nameCi,runTime.timeName(),mesh,IOobject::READ_IF_PRESENT,IOobject::AUTO_WRITE),
			 mesh,
			 dimensionedScalar("",dimMoles/dimVol,0.), //freak.solu_conc[i]
			 bcTypes
			)
		);
}

Info <<"end ECwifields " << endl;