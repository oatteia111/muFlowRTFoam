dimensionedScalar cpw(transportProperties.lookupOrDefault("cpw",dimensionedScalar("cpw",dimensionSet(0,2,-2,-1,0,0,0),4182)));
dimensionedScalar lbdaTw(transportProperties.lookupOrDefault("lbdaTw",dimensionedScalar("lbdaTw",dimensionSet(1,1,-3,-1,0,0,0),0.6)));
scalar activateThermal(transportProperties.lookupOrDefault<scalar>("activateThermal",0));
scalar dT=1e-9;scalar dT1=1e-9;

//NB (dimensionType MASS , LENGTH , TIME , TEMPERATURE ,MOLES , CURRENT , LUMINOUS_INTENSITY)
Info<< "Reading field T\n" << endl;
volScalarField T  // Temperature
(
    IOobject("T",runTime.timeName(),mesh,IOobject::READ_IF_PRESENT,IOobject::AUTO_WRITE),mesh,dimensionedScalar("T",dimTemperature,25.)
);


volScalarField cps  // heat capacity of the solid [L2 T-2 °C-1]
(
	IOobject("cps",runTime.constant(),mesh,IOobject::READ_IF_PRESENT,IOobject::AUTO_WRITE),mesh,dimensionedScalar("cps",dimensionSet(0,2,-2,-1,0,0,0),840)
);

volScalarField lbdaTs  // Thermal conductivity of the solid [M L T-3 °C-1]
(
	IOobject("lbdaTs",runTime.constant(),mesh,IOobject::READ_IF_PRESENT,IOobject::AUTO_WRITE),mesh,dimensionedScalar("lbdaTs",dimensionSet(1,1,-3,-1,0,0,0),3.5)
);

volScalarField lbdaT = sw*eps*lbdaTw + (1-eps)*lbdaTs;


// creating the internal Su and Sp fields to input temperature conditions inside the domain
volScalarField::Internal TSu
(
   IOobject("TSu",runTime.timeName(),mesh,IOobject::NO_READ,IOobject::NO_WRITE),
		mesh,
		dimensionedScalar("zero",dimensionSet(0,0,-1,1,0,0,0),0.0) // the overall eq in in K.s-1
		//dimensioned<scalar>("zero",eqn.dimensions()/dimVolume/psi.dimensions(),0.0), false
);

volScalarField::Internal TSp
(
   IOobject("TSp",runTime.timeName(),mesh,IOobject::NO_READ,IOobject::NO_WRITE),
		mesh,
		dimensionedScalar("zero",dimensionSet(0,0,-1,0,0,0,0),0.0) 
		//dimensioned<scalar>("zero",eqn.dimensions()/dimVolume/psi.dimensions(),0.0), false
);

labelList cellsTfix0; std::vector<float> condTfix;
labelList cellsTwel0; std::vector<float> condTwel;

if (fexists(cur_dir+"/constant/options/tfix")) {cellSet selectedTfix(mesh, "tfix");cellsTfix0 = selectedTfix.toc();}
if (fexists(cur_dir+"/constant/options/twel")) {cellSet selectedTwel(mesh, "twel");cellsTwel0 = selectedTwel.toc();}

int itTfix=0;output = getCbuffer("tfix",itTfix,cellsTfix0.size()); 
float tTfix = output.t; std::vector<float> cellsTfixData=output.d;

int itTwel=0;output = getCbuffer("twel",itTwel,cellsTwel0.size()); 
float tTwel = output.t; std::vector<float> cellsTwelData=output.d;

condTfix.resize(cellsTfix0.size());
condTwel.resize(cellsTwel0.size());
std::vector<int> cellsTfix = indexC(cellsTfix0,cellsTfixData) ;
std::vector<int> cellsTwel = indexC(cellsTwel0,cellsTwelData) ;
Info<< "All T read\n" << endl;
