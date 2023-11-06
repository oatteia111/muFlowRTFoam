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
labelList cellsCfix0; std::vector<float> condCfix;std::vector<float> cellsCfixData;
labelList cellsCwel0; std::vector<float> condCwel;std::vector<float> cellsCwelData;
labelList cellsCghb0; std::vector<float> condCghb;std::vector<float> cellsCghbData;
labelList cellsCrch0; std::vector<float> condCrch;std::vector<float> cellsCrchData;
int itCfix=0;int itCwel=0;int itCghb=0;int itCrch=0;
float tCwel=0;float tCfix=0;float tCghb=0; float tCrch=0;

if (fexists(cur_dir+"/constant/options/cfix")) {cellSet selectedCfix(mesh, "cfix");cellsCfix0 = selectedCfix.toc();}
if (fexists(cur_dir+"/constant/options/cwel")) {cellSet selectedCwel(mesh, "cwel");cellsCwel0 = selectedCwel.toc();}
if (fexists(cur_dir+"/constant/options/cghb")) {cellSet selectedCghb(mesh, "cghb");cellsCghb0 = selectedCghb.toc();}
if (fexists(cur_dir+"/constant/options/crch")) {cellSet selectedCrch(mesh, "crch");cellsCrch0 = selectedCrch.toc();}

while (time>=tCfix) {output = getCbuffer("cfix",itCfix,cellsCfix0.size());itCfix+=1;tCfix = output.t;cellsCfixData=output.d;}
while (time>=tCwel) {output = getCbuffer("cwel",itCwel,cellsCwel0.size());itCwel+=1;tCwel = output.t;cellsCwelData=output.d;}
while (time>=tCghb) {output = getCbuffer("cghb",itCghb,cellsCghb0.size());itCghb+=1;tCghb = output.t;cellsCghbData=output.d;}
while (time>=tCrch) {output = getCbuffer("crch",itCrch,cellsCrch0.size());itCrch+=1;tCrch = output.t;cellsCrchData=output.d;}

condCfix.resize(cellsCfix0.size());
condCwel.resize(cellsCwel0.size());
condCghb.resize(cellsCghb0.size());
condCrch.resize(cellsCrch0.size());

std::vector<int> cellsCfix = indexC(cellsCfix0,cellsCfixData) ;
std::vector<int> cellsCwel = indexC(cellsCwel0,cellsCwelData) ;
std::vector<int> cellsCghb = indexC(cellsCghb0,cellsCghbData) ;
std::vector<int> cellsCrch = indexC(cellsCrch0,cellsCrchData) ;

//for (i=0;i<20;i++) {Info<<"rch "<< i<<" conc "<<cellsCrchData[i*4+2]<<endl;}

