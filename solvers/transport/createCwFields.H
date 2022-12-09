dimensionedScalar alphaL(transportProperties.lookupOrDefault("alphaL",dimensionedScalar("alphaL",dimLength,0.0)));
dimensionedScalar alphaT(transportProperties.lookupOrDefault("alphaT",dimensionedScalar("alphaT",dimLength,0.0)));
dimensionedScalar lbdaw(transportProperties.lookupOrDefault("lbdaw",dimensionedScalar("lbdaw",dimless/dimTime,0.0)));
dimensionedScalar rhos(transportProperties.lookupOrDefault("rhos",dimensionedScalar("rhos",dimless,2.65)));
dimensionedScalar Dw0(transportProperties.lookupOrDefault("Dw0",dimensionedScalar("Dw0",dimArea/dimTime,1e-10)));
scalar transportSteady(transportProperties.lookupOrDefault<scalar>("transportSteady",0));
//const scalar dCmax = runTime.controlDict().lookupOrDefault("dCmax", 0.01);
//const scalar dCresidual = runTime.controlDict().lookupOrDefault("dCresidual", 1e-3);
Info<< "Reading field Cw\n" << endl;
wordList bcTypes (Kh.boundaryField().size(),zeroGradientFvPatchScalarField::typeName);

volScalarField Cw
(
    IOobject("Cw",runTime.timeName(),mesh,IOobject::READ_IF_PRESENT,IOobject::AUTO_WRITE),mesh,dimensionedScalar("",dimMass/dimVol,0),bcTypes
);

volScalarField Kd
(
    IOobject("Kd",runTime.constant(),mesh,IOobject::READ_IF_PRESENT,IOobject::AUTO_WRITE),mesh,dimensionedScalar("",dimless,0),bcTypes
);

volTensorField Deffw
(
    IOobject("Deffw",runTime.timeName(),mesh,IOobject::READ_IF_PRESENT,IOobject::NO_WRITE),mesh,dimensionedTensor("",dimArea/dimTime,Zero) 
);

volScalarField R = 1 + (1-eps)/eps*rhos*Kd;

Info<< "All C read\n" << endl;
