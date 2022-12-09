const dimensionedScalar  epsilonK ("vacuumPermittivity",dimensionSet(-1, -3, 4, 0, 0, 2, 0),8.8541878176e-12  ); // ok valid for units and value
const dimensionedScalar AK("AvogradoNumber", dimensionSet(0, 0, 0, 0, -1, 0, 0), 6.022140857e+23);
const dimensionedScalar eK("elementaryCharge", dimensionSet( 0, 0, 1, 0, 0, 1, 0 ),1.6021766208e-19); // in coulomb = A.s
const dimensionedScalar kbK( "BoltzmannConstant", dimensionSet( 1, 2, -2, -1, 0, 0, 0 ),1.38064852e-23);
const dimensionedScalar FK( "FaradayConstant", (eK*AK).dimensions(),(eK*AK).value() ); // ok 99680 C.mol-1 or A.s.mol-1 [0,0,1,0-1,1,0] 

dimensionedScalar D25("diffCoeff", dimensionSet(0, 2, -1, 0, 0, 0, 0), 1e-12);
dimensionedVector extraE (transportProperties.lookupOrDefault<dimensionedVector>("extraEField", dimensionedVector("0", dimensionSet(1, 1, -3, 0, 0, -1, 0), vector::zero)));
scalar relPerm (transportProperties.lookupOrDefault<scalar>( "dielectricConstant", 1. ));
dimensionedScalar TK(transportProperties.lookupOrDefault("tempKelvin",dimensionedScalar("tempKelvin",dimTemperature,298)));
volScalarField psiE
(
    IOobject("psiE",runTime.timeName(),mesh,IOobject::READ_IF_PRESENT,IOobject::AUTO_WRITE),
	mesh,
	dimensionedScalar("psiE",dimensionSet(1,2,-3,0,0,-1,0),0)  //this is in volts=kg.m^2/(A.s^3)
);
volScalarField psiI
(
    IOobject("psiI",runTime.timeName(),mesh,IOobject::READ_IF_PRESENT,IOobject::AUTO_WRITE),
	mesh,
	dimensionedScalar("psiI",dimensionSet(1,2,-3,0,0,-1,0),0)
);
