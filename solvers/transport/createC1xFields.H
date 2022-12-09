Info<< "Reading C1x transport props\n" << endl;

dimensionedScalar Dg0(transportProperties.lookupOrDefault("Dg0",dimensionedScalar("Dg0",dimArea/dimTime,1e-6)));
scalar Kdg(transportProperties.lookupOrDefault<scalar>("Kdg",0.));
scalar Kds(transportProperties.lookupOrDefault<scalar>("Kds",0.));
dimensionedScalar kmg(transportProperties.lookupOrDefault("kmg",dimensionedScalar("kmg",dimless/dimTime,0.0)));
dimensionedScalar kms(transportProperties.lookupOrDefault("kms",dimensionedScalar("kms",dimless/dimTime,0.0)));
dimensionedScalar lbdaw(transportProperties.lookupOrDefault("lbdaw",dimensionedScalar("lbdaw",dimless/dimTime,0.0)));
dimensionedScalar lbdag(transportProperties.lookupOrDefault("lbdag",dimensionedScalar("lbdag",dimless/dimTime,0.0)));
dimensionedScalar lbdas(transportProperties.lookupOrDefault("lbdas",dimensionedScalar("lbdas",dimless/dimTime,0.0)));

volScalarField Cw
(
    IOobject("Cw",runTime.timeName(), mesh,IOobject::READ_IF_PRESENT,IOobject::AUTO_WRITE),mesh,dimensionedScalar("",dimMass/dimVol,0) 
);
volScalarField Cg
(
    IOobject("Cg",runTime.timeName(), mesh,IOobject::READ_IF_PRESENT,IOobject::AUTO_WRITE),mesh,dimensionedScalar("",dimMass/dimVol,0) 
);
volScalarField Cs
(
    IOobject("Cs",runTime.timeName(),mesh,IOobject::READ_IF_PRESENT,IOobject::NO_WRITE),mesh,dimensionedScalar("",dimMass/dimVol,0) 
);

volTensorField Deffg
(
    IOobject("Deffg",runTime.timeName(),mesh,IOobject::READ_IF_PRESENT,IOobject::NO_WRITE),mesh,dimensionedTensor("",dimArea/dimTime,Zero) 
);
dimensionedScalar Deffw(transportProperties.lookupOrDefault("Deffw",dimensionedScalar("Deffw",dimArea/dimTime,0)));

volScalarField Defft ("Defft", unity*1e-10);
volScalarField sg ("sg", 1-sw);
volScalarField thet1 ("thet1", unity);


volScalarField massTransfer_g
(
    IOobject("massTransfer_g",runTime.timeName(),mesh,IOobject::READ_IF_PRESENT,IOobject::NO_WRITE),mesh,dimensionedScalar("",dimMass/(dimVol*dimTime),0.) 
);
volScalarField massTransfer_s
(
    IOobject("massTransfer_s",runTime.timeName(),mesh,IOobject::READ_IF_PRESENT,IOobject::NO_WRITE),mesh,dimensionedScalar("",dimMass/(dimVol*dimTime),0.) 
);
