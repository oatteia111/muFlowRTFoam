dimensionedScalar alphaL(transportProperties.lookupOrDefault("alphaL",dimensionedScalar("alphaL",dimLength,0.0)));
dimensionedScalar alphaT(transportProperties.lookupOrDefault("alphaT",dimensionedScalar("alphaT",dimLength,0.0)));
dimensionedScalar lbdaw(transportProperties.lookupOrDefault("lbdaw",dimensionedScalar("lbdaw",dimless/dimTime,0.0)));
dimensionedScalar Dw0(transportProperties.lookupOrDefault("Dw0",dimensionedScalar("Dw0",dimArea/dimTime,1e-10)));
dimensionedScalar Koc(transportProperties.lookupOrDefault("Koc",dimensionedScalar("Koc",dimVol/dimMass,0)));
scalar transportSteady(transportProperties.lookupOrDefault<scalar>("transportSteady",0));
double dC,mxC,dff,dC0,sC,sC1,dtForC,dtForChem;
const scalar dCmax = runTime.controlDict().lookupOrDefault("dCmax", 0.02);
const scalar dTmax = runTime.controlDict().lookupOrDefault("dTmax", 0.02);

//const scalar dCresidual = runTime.controlDict().lookupOrDefault("dCresidual", 1e-3);
Info<< "Reading field C\n" << endl;
wordList bcTypes (Kh.boundaryField().size(),zeroGradientFvPatchScalarField::typeName);
volScalarField rhos
(
    IOobject("rhos",runTime.constant(),mesh,IOobject::READ_IF_PRESENT,IOobject::NO_WRITE),mesh,dimensionedScalar("rhos",dimMass/dimVol,2650),bcTypes
);
volScalarField C
(
    IOobject("C",runTime.timeName(),mesh,IOobject::READ_IF_PRESENT,IOobject::AUTO_WRITE),mesh,dimensionedScalar("C",dimMass/dimVol,0),bcTypes
);

volScalarField foc //organic matter content
(
	IOobject("foc",runTime.constant(),mesh,IOobject::READ_IF_PRESENT,IOobject::AUTO_WRITE),mesh,dimensionedScalar("foc",dimless,0),bcTypes
);
volScalarField dgrain //grain diameter (mm), used for PFAs sorption
(IOobject("dgrain",runTime.constant(),mesh,IOobject::READ_IF_PRESENT,IOobject::NO_WRITE),mesh,dimensionedScalar("dgrain",dimless,1e-2),bcTypes);
volScalarField Rf //retardation factor (can be transient for AW sorption)
(
	IOobject("Rf",runTime.timeName(),mesh,IOobject::READ_IF_PRESENT,IOobject::NO_WRITE),mesh,dimensionedScalar("Rf",dimless,1),bcTypes
);
Rf = 1 + (1-eps)/eps*rhos*Koc*foc;
//scalar activateSorptionAW(transportProperties.lookupOrDefault<scalar>("activateSorptionAW",0));std::cout<<"AW "<<activateSorptionAW<<"\n";
volTensorField Deffw
(
    IOobject("Deffw",runTime.timeName(),mesh,IOobject::READ_IF_PRESENT,IOobject::NO_WRITE),mesh,dimensionedTensor("Deffw",dimArea/dimTime,Zero) 
);

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
labelList cellsCfix0; std::vector<float> cellsCfixData;float tCfix=0;int itCfix=0;
labelList cellsCwel0; std::vector<float> cellsCwelData;float tCwel=0;int itCwel=0;
labelList cellsCriv0; std::vector<float> cellsCrivData;float tCriv=0;int itCriv=0; 
labelList cellsCghb0; std::vector<float> cellsCghbData;float tCghb=0;int itCghb=0; 
labelList cellsCrch0; std::vector<float> cellsCrchData;float tCrch=0;int itCrch=0;

// reading the data for hfix
//std::ifstream inputCfxdata{cur_dir+"/constant/options/hfix"}; //
//cellsCfxData = {std::istream_iterator<double>{inputCfxdata}, {}};
if (fexists(cur_dir+"/constant/options/cfix")) { cellSet selectedCfix(mesh, "cfix");cellsCfix0 = selectedCfix.toc();}
if (fexists(cur_dir+"/constant/options/cwel")) { cellSet selectedCwel(mesh, "cwel");cellsCwel0 = selectedCwel.toc();}
if (fexists(cur_dir+"/constant/options/cghb")) { cellSet selectedCghb(mesh, "cghb");cellsCghb0 = selectedCghb.toc();}
if (fexists(cur_dir+"/constant/options/criv")) { cellSet selectedCriv(mesh, "criv");cellsCriv0 = selectedCriv.toc();}
if (fexists(cur_dir+"/constant/options/crch")) { cellSet selectedCrch(mesh, "crch");cellsCrch0 = selectedCrch.toc();}
//get the data
time = mesh.time().value()/day;
while (time>=tCfix) {output = getCbuffer("cfix",itCfix,cellsCfix0.size());itCfix +=1;tCfix = output.t;cellsCfixData=output.d;tnext=tCfix;}
while (time>=tCwel) {output = getCbuffer("cwel",itCwel,cellsCwel0.size());itCwel +=1;tCwel = output.t;cellsCwelData=output.d;tnext=tCwel;}
while (time>=tCghb) {output = getCbuffer("cghb",itCghb,cellsCghb0.size());itCghb +=1;tCghb = output.t;cellsCghbData=output.d;tnext=tCghb;}
while (time>=tCriv) {output = getCbuffer("criv",itCriv,cellsCriv0.size());itCriv +=1;tCriv = output.t;cellsCrivData=output.d;tnext=tCriv;}
while (time>=tCrch) {output = getCbuffer("crch",itCrch,cellsCrch0.size());itCrch +=1;tCrch = output.t;cellsCrchData=output.d;tnext=tCrch;}

std::vector<int> cellsCfix = indexC(cellsCfix0,cellsCfixData) ;
std::vector<int> cellsCwel = indexC(cellsCwel0,cellsCwelData) ;
std::vector<int> cellsCghb = indexC(cellsCghb0,cellsCghbData) ;
std::vector<int> cellsCriv = indexC(cellsCriv0,cellsCrivData) ;
std::vector<int> cellsCrch = indexC(cellsCrch0,cellsCrchData) ;

