dimensionedScalar alphaL(transportProperties.lookupOrDefault("alphaL",dimensionedScalar("alphaL",dimLength,0.0)));
dimensionedScalar alphaT(transportProperties.lookupOrDefault("alphaT",dimensionedScalar("alphaT",dimLength,0.0)));
dimensionedScalar lbdaw(transportProperties.lookupOrDefault("lbdaw",dimensionedScalar("lbdaw",dimless/dimTime,0.0)));
dimensionedScalar Dw0(transportProperties.lookupOrDefault("Dw0",dimensionedScalar("Dw0",dimArea/dimTime,1e-10)));
scalar transportSteady(transportProperties.lookupOrDefault<scalar>("transportSteady",0));
//const scalar dCmax = runTime.controlDict().lookupOrDefault("dCmax", 0.01);
//const scalar dCresidual = runTime.controlDict().lookupOrDefault("dCresidual", 1e-3);
Info<< "Reading field C\n" << endl;
wordList bcTypes (Kh.boundaryField().size(),zeroGradientFvPatchScalarField::typeName);

volScalarField C
(
    IOobject("C",runTime.timeName(),mesh,IOobject::READ_IF_PRESENT,IOobject::AUTO_WRITE),mesh,dimensionedScalar("C",dimMass/dimVol,0),bcTypes
);

volScalarField rhos
(
    IOobject("rhos",runTime.constant(),mesh,IOobject::READ_IF_PRESENT,IOobject::AUTO_WRITE),mesh,dimensionedScalar("rhos",dimMass/dimVol,0),bcTypes
);

volScalarField Kds //it was impossible to call it Kd no dimension.. why???
(
    IOobject("Kds",runTime.constant(),mesh,IOobject::READ_IF_PRESENT,IOobject::AUTO_WRITE),mesh,dimensionedScalar("Kds",dimensionSet(-1,3,0,0,0,0,0),0),bcTypes
);

volTensorField Deffw
(
    IOobject("Deffw",runTime.timeName(),mesh,IOobject::READ_IF_PRESENT,IOobject::NO_WRITE),mesh,dimensionedTensor("Deffw",dimArea/dimTime,Zero) 
);
//Info<< " dim rhos"<<rhos.dimensions()<<" Kds "<<Kds.dimensions()<< endl;
volScalarField R = 1 + (1-eps)/eps*rhos*Kds;

Info<< "All C read\n" << endl;
