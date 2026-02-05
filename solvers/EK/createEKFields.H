scalar activateEK(transportProperties.lookupOrDefault<scalar>("activateEK",0));

const dimensionedScalar  epsilonK ("vacuumPermittivity",dimensionSet(-1, -3, 4, 0, 0, 2, 0),8.8541878176e-12  ); // ok valid for units and value
const dimensionedScalar AK("AvogradoNumber", dimensionSet(0, 0, 0, 0, -1, 0, 0), 6.022140857e+23);
const dimensionedScalar eK("elementaryCharge", dimensionSet( 0, 0, 1, 0, 0, 1, 0 ),1.6021766208e-19); // in coulomb = A.s
//const dimensionedScalar kbK( "BoltzmannConstant", dimensionSet( 1, 2, -2, -1, 0, 0, 0 ),1.38064852e-23);
const dimensionedScalar FK( "FaradayConstant", (eK*AK).dimensions(),(eK*AK).value() ); // ok 99680 C.mol-1 or A.s.mol-1 [0,0,1,0-1,1,0] 
const dimensionedScalar  Rgas ("gasConstant",dimensionSet(1,2,-2,-1,-1,0,0),8.31446  ); // ok valid for units and value

dimensionedScalar D25("D25", dimensionSet(0,2,-1,0,0,0,0), 1e-12);
//volScalarField D25(IOobject("D25",runTime.timeName(),mesh,IOobject::READ_IF_PRESENT,IOobject::AUTO_WRITE),mesh,dimensionedScalar("D25", dimensionSet(0,2,-1,0,0,0,0), 1e-12));
dimensionedVector extraE (transportProperties.lookupOrDefault<dimensionedVector>("extraEField", dimensionedVector("0", dimensionSet(1, 1, -3, 0, 0, -1, 0), vector::zero)));
scalar relPerm (transportProperties.lookupOrDefault<scalar>( "dielectricConstant", 80. )); //water relative permittivity (seems seawater is 81)
dimensionedScalar TK(transportProperties.lookupOrDefault("tempKelvin",dimensionedScalar("tempKelvin",dimTemperature,298)));

volScalarField JE(IOobject("JE",runTime.timeName(),mesh,IOobject::READ_IF_PRESENT,IOobject::AUTO_WRITE),mesh,dimensionedScalar("JE",dimensionSet(1,-2,-1,0,0,0,0),0.));
volScalarField JEpos = JE;
volScalarField JEneg = JE;
volScalarField fJE(IOobject("fJE",runTime.timeName(),mesh,IOobject::READ_IF_PRESENT,IOobject::AUTO_WRITE),mesh,dimensionedScalar("fJE",dimless,0.));
volScalarField fJE1 = fJE;
volScalarField JE1 = JE;
volScalarField JE_tmp = JE;
volScalarField cf = D25 * FK / (Rgas*TK)*eps;
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
	dimensionedScalar("psiI",dimensionSet(1,2,-3,0,0,-1,0),0),
	bcTypes
);

dimensionedScalar cum("0", dimensionSet(0, -3, 0, 0, 1, 0, 0), 0);   
volScalarField rhoE( psiE * dimensionedScalar("norm", dimless/dimArea, 0.) );
volScalarField souE = rhoE*1;         
volScalarField souI = rhoE*1;
dimensionedScalar invMolweight("iweigth",dimensionSet(-1, 0, 0, 0, 1, 0, 0),1);

// -------------- conditions for EK -------------
labelList cellsEwel0; std::vector<float> qEwel;std::vector<float> cellsEwelData;float tEwel=0;int itEwel=0;
if (fexists(cur_dir+"/constant/options/ewel")) { 
cellSet selectedEwel(mesh, "ewel");cellsEwel0 = selectedEwel.toc();
while (time>=tEwel) {output = getCbuffer("ewel",itEwel,cellsEwel0.size());itEwel +=1;tEwel = output.t;tnext=min(tnext,tEwel);cellsEwelData=output.d;}
}