dimensionedScalar alphaL(transportProperties.lookupOrDefault("alphaL",dimensionedScalar("alphaL",dimLength,0.0)));
dimensionedScalar alphaT(transportProperties.lookupOrDefault("alphaT",dimensionedScalar("alphaT",dimLength,0.0)));
dimensionedScalar lbdaw(transportProperties.lookupOrDefault("lbdaw",dimensionedScalar("lbdaw",dimless/dimTime,0.0)));
dimensionedScalar Dw0(transportProperties.lookupOrDefault("Dw0",dimensionedScalar("Dw0",dimArea/dimTime,1e-10)));
scalar transportSteady(transportProperties.lookupOrDefault<scalar>("transportSteady",0));
double dC=0;double mxC=0;double sC=0;double sC1=0;double dtForC;
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
    IOobject("rhos",runTime.constant(),mesh,IOobject::READ_IF_PRESENT,IOobject::AUTO_WRITE),mesh,dimensionedScalar("rhos",dimMass/dimVol,2650),bcTypes
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

// creating the internal Su and Sp fields to input temperature conditions inside the domain
volScalarField::Internal CSu
(
   IOobject("CSu",runTime.timeName(),mesh,IOobject::NO_READ,IOobject::NO_WRITE),
		mesh,
		dimensionedScalar("zero",dimensionSet(1,-3,-1,0,0,0,0),0.0) 
);

volScalarField::Internal CSp
(
   IOobject("CSp",runTime.timeName(),mesh,IOobject::NO_READ,IOobject::NO_WRITE),
		mesh,
		dimensionedScalar("zero",dimensionSet(0,0,-1,0,0,0,0),0.0) 
);

// -------------- conditions for concentration -------------
labelList cellsCfix0; std::vector<float> cellsCfixData;std::vector<float> condCfix;
labelList cellsCwel0; std::vector<float> cellsCwelData;std::vector<float> condCwel;
labelList cellsCghb0; std::vector<float> cellsCghbData;std::vector<float> condCghb;
labelList cellsCrch0; std::vector<float> cellsCrchData;std::vector<float> condCrch;

// reading the data for hfix
std::ifstream inputCfixdata{cur_dir+"/constant/options/cfix"}; //
cellsCfixData = {std::istream_iterator<float>{inputCfixdata}, {}};
if (cellsCfixData.size()>0) { cellSet selectedCfix(mesh, "cfix"); cellsCfix0 = selectedCfix.toc();}
condCfix.resize(cellsCfix0.size());
std::vector<int> cellsCfix = indexC(cellsCfix0,cellsCfixData) ;

// reading the cwell data for wells
std::ifstream inputCweldata{cur_dir+"/constant/options/cwel"}; //
cellsCwelData = {std::istream_iterator<float>{inputCweldata}, {}};
if (cellsCwelData.size()>0) {	cellSet selectedCwel(mesh, "cwel");cellsCwel0 = selectedCwel.toc();}
condCwel.resize(cellsCwel0.size());
std::vector<int> cellsCwel = indexC(cellsCwel0,cellsCwelData) ;
//for (i=0;i<cellsCw.size();i++) {Info<<"i "<<cellsCw[i]<<" center "<<mesh.C()[cellsCw[i]]<<endl;}

// reading the data for ghb
std::ifstream inputCghbdata{cur_dir+"/constant/options/cghb"}; //
cellsCghbData = {std::istream_iterator<float>{inputCghbdata}, {}};
if (cellsCghbData.size()>0) { cellSet selectedCghb(mesh, "cghb"); cellsCghb0 = selectedCghb.toc();}
condCghb.resize(cellsCghb0.size());
std::vector<int> cellsCghb = indexC(cellsCghb0,cellsCghbData) ;

// reading the data for rch
std::ifstream inputCrchdata{cur_dir+"/constant/options/crch"}; //
cellsCrchData = {std::istream_iterator<float>{inputCrchdata}, {}};
if (cellsCrchData.size()>0) { cellSet selectedCrch(mesh, "crch"); cellsCrch0 = selectedCrch.toc();}
condCrch.resize(cellsCrch0.size());
std::vector<int> cellsCrch = indexC(cellsCrch0,cellsCrchData) ;

