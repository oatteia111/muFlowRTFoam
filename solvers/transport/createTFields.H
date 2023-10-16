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

//////////////////////// read conditions //////////////////////////
labelList cellsTfix0; std::vector<float> condTfix;std::vector<float> cellsTfixData;
labelList cellsTwel0; std::vector<float> condTwel;std::vector<float> cellsTwelData;
labelList cellsTghb0; std::vector<float> condTghb;std::vector<float> cellsTghbData;
labelList cellsTrch0; std::vector<float> condTrch;std::vector<float> cellsTrchData;
int itTfix=0;int itTwel=0;int itTghb=0;int itTrch=0;
float tTwel=0;float tTfix=0;float tTghb=0; float tTrch=0;

if (fexists(cur_dir+"/constant/options/tfix")) {cellSet selectedTfix(mesh, "tfix");cellsTfix0 = selectedTfix.toc();}
if (fexists(cur_dir+"/constant/options/twel")) {cellSet selectedTwel(mesh, "twel");cellsTwel0 = selectedTwel.toc();}
if (fexists(cur_dir+"/constant/options/tghb")) {cellSet selectedTghb(mesh, "tghb");cellsTghb0 = selectedTghb.toc();}
if (fexists(cur_dir+"/constant/options/trch")) {cellSet selectedTrch(mesh, "trch");cellsTrch0 = selectedTrch.toc();}

while (time>=tTfix) {output = getCbuffer("tfix",itTfix,cellsTfix0.size());itTfix+=1;tTfix = output.t;cellsTfixData=output.d;}
while (time>=tTwel) {output = getCbuffer("twel",itTwel,cellsTwel0.size());itTwel+=1;tTwel = output.t;cellsTwelData=output.d;}
while (time>=tTghb) {output = getCbuffer("tghb",itTghb,cellsTghb0.size());itTghb+=1;tTghb = output.t;cellsTghbData=output.d;}
while (time>=tTrch) {output = getCbuffer("trch",itTrch,cellsTrch0.size());itTrch+=1;tTrch = output.t;cellsTrchData=output.d;}

condTfix.resize(cellsTfix0.size());
condTwel.resize(cellsTwel0.size());
condTghb.resize(cellsTghb0.size());
condTrch.resize(cellsTrch0.size());

std::vector<int> cellsTfix = indexC(cellsTfix0,cellsTfixData) ;
std::vector<int> cellsTwel = indexC(cellsTwel0,cellsTwelData) ;
std::vector<int> cellsTghb = indexC(cellsTghb0,cellsTghbData) ;
std::vector<int> cellsTrch = indexC(cellsTrch0,cellsTrchData) ;
Info<< "All T read\n" << endl;
