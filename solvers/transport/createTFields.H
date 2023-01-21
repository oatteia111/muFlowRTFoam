dimensionedScalar cpw(transportProperties.lookupOrDefault("cpw",dimensionedScalar("cpw",dimensionSet(0,2,-2,-1,0,0,0),4182)));
dimensionedScalar lbdaTw(transportProperties.lookupOrDefault("lbdaTw",dimensionedScalar("lbdaTw",dimensionSet(1,1,-3,-1,0,0,0),0.6)));
scalar activateThermal(transportProperties.lookupOrDefault<scalar>("activateThermal",0));

//NB (dimensionType MASS , LENGTH , TIME , TEMPERATURE ,MOLES , CURRENT , LUMINOUS_INTENSITY)
Info<< "Reading field T\n" << endl;
volScalarField T  // Temperature
(
    IOobject("T",runTime.timeName(),mesh,IOobject::READ_IF_PRESENT,IOobject::AUTO_WRITE),mesh,dimensionedScalar("T",dimTemperature,25.)
);
volScalarField cps  // heat capacity of the solid [L2 T-2 Θ-1]
(
    IOobject("cps",runTime.constant(),mesh,IOobject::READ_IF_PRESENT,IOobject::AUTO_WRITE),mesh,dimensionedScalar("cps",dimensionSet(0,2,-2,-1,0,0,0),840)
);

volScalarField lbdaTs  // Thermal conductivity of the solid [M L T-3 Θ-1]
(
    IOobject("lbdaTs",runTime.constant(),mesh,IOobject::READ_IF_PRESENT,IOobject::AUTO_WRITE),mesh,dimensionedScalar("lbdaTs",dimensionSet(1,1,-3,-1,0,0,0),3.5)
);

volScalarField lbdaT = sw*eps*lbdaTw + (1-eps)*lbdaTs;

//options
// reading the hwell data required for injection of temperature (cells Hw and tw shall be the same)
cellSet selectedCells(mesh, "hwel");
labelList cellsHw = selectedCells.toc();
std::ifstream inputHwdata{cur_dir+"/constant/options/hwel"}; //
std::vector<float> cellsHwData;
cellsHwData = {std::istream_iterator<float>{inputHwdata}, {}};
//reading the twell data
cellSet selectedCells1(mesh, "twel");
labelList cellsTw = selectedCells1.toc();
std::ifstream inputTwdata{cur_dir+"/constant/options/twel"}; //
std::vector<float> cellsTwData;
cellsTwData = {std::istream_iterator<float>{inputTwdata}, {}};

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


Info<< "All T read\n" << endl;
