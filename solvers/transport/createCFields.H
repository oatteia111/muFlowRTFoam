dimensionedScalar alphaL(transportProperties.lookupOrDefault("alphaL",dimensionedScalar("alphaL",dimLength,0.0)));
dimensionedScalar alphaT(transportProperties.lookupOrDefault("alphaT",dimensionedScalar("alphaT",dimLength,0.0)));
dimensionedScalar lbdaw(transportProperties.lookupOrDefault("lbdaw",dimensionedScalar("lbdaw",dimless/dimTime,0.0)));
dimensionedScalar Dw0(transportProperties.lookupOrDefault("Dw0",dimensionedScalar("Dw0",dimArea/dimTime,1e-10)));
dimensionedScalar Koc(transportProperties.lookupOrDefault("Koc",dimensionedScalar("Koc",dimVol/dimMass,0)));
scalar transportSteady(transportProperties.lookupOrDefault<scalar>("transportSteady",0));
double dC,mnC,mxC,dff,dC0,sC,sC1;
const scalar dCmax = runTime.controlDict().lookupOrDefault("dCmax", 0.02);
const scalar dTmax = runTime.controlDict().lookupOrDefault("dTmax", 0.02);

scalar reactIndicator = 0; // indicator for reaction; 0 crank; 1 euler
//adaptiveReactiveDdtScheme<scalar> Rddt(mesh, "Cwi", phaseIndicator);

//const scalar dCresidual = runTime.controlDict().lookupOrDefault("dCresidual", 1e-3);
Info<< "Reading field C\n" << endl;
wordList bcTypes (Kh.boundaryField().size(),zeroGradientFvPatchScalarField::typeName);
volScalarField rhos(IOobject("rhos",runTime.constant(),mesh,IOobject::READ_IF_PRESENT,IOobject::NO_WRITE),mesh,dimensionedScalar("rhos",dimMass/dimVol,2650),bcTypes);
volScalarField C(IOobject("C",runTime.timeName(),mesh,IOobject::READ_IF_PRESENT,IOobject::AUTO_WRITE),mesh,dimensionedScalar("C",dimMass/dimVol,0),bcTypes);

volScalarField foc //organic matter content
(
	IOobject("foc",runTime.constant(),mesh,IOobject::READ_IF_PRESENT,IOobject::AUTO_WRITE),mesh,dimensionedScalar("foc",dimless,0),bcTypes
);
//volScalarField dgrain //grain diameter (mm), used for PFAs sorption
//(IOobject("dgrain",runTime.constant(),mesh,IOobject::READ_IF_PRESENT,IOobject::NO_WRITE),mesh,dimensionedScalar("dgrain",dimless,1e-2),bcTypes);
volScalarField Rf
(	IOobject("Rf",runTime.timeName(),mesh,IOobject::READ_IF_PRESENT,IOobject::NO_WRITE),mesh,dimensionedScalar("Rf",dimless,1),bcTypes);
//Rf=eps*2;Info<< "Rf "<<Rf[0] <<"\n" << endl;
//Rf = 1 + (1-eps)/eps*rhos*Koc*foc;
//scalar activateSorptionAW(transportProperties.lookupOrDefault<scalar>("activateSorptionAW",0));std::cout<<"AW "<<activateSorptionAW<<"\n";
volTensorField Deffw(IOobject("Deffw",runTime.timeName(),mesh,IOobject::READ_IF_PRESENT,IOobject::NO_WRITE),mesh,dimensionedTensor("Deffw",dimArea/dimTime,Zero) );

// creating the internal Su and Sp fields to input temperature conditions inside the domain
volScalarField::Internal CSu
(
   IOobject("CSu",runTime.timeName(),mesh,IOobject::NO_READ,IOobject::NO_WRITE),
		mesh,dimensionedScalar("zero",dimensionSet(1,-3,-1,0,0,0,0),0.0) 
);

volScalarField::Internal CSp
(
   IOobject("CSp",runTime.timeName(),mesh,IOobject::NO_READ,IOobject::NO_WRITE),
		mesh,dimensionedScalar("zero",dimensionSet(0,0,-1,0,0,0,0),0.0) 
);

// -------------- conditions for concentration -------------
// thre is chfix: conc for the cells where h is fixed and cfix: ccells for which only c is fixed
labelList cellsChfix0; std::vector<float> cellsChfixData;float tChfix=0;int itChfix=0;
labelList cellsCfix0; std::vector<float> cellsCfixData;float tCfix=0;int itCfix=0;
labelList cellsCwel0; std::vector<float> cellsCwelData;float tCwel=0;int itCwel=0;
labelList cellsCriv0; std::vector<float> cellsCrivData;float tCriv=0;int itCriv=0; 
labelList cellsCghb0; std::vector<float> cellsCghbData;float tCghb=0;int itCghb=0; 
labelList cellsCrch0; std::vector<float> cellsCrchData;float tCrch=0;int itCrch=0;

// reading the data for hfix
time = mesh.time().value() ;
std::vector<int> cellsChfix,cellsCfix,cellsCwel,cellsCghb,cellsCriv,cellsCrch;

if (fexists(cur_dir+"/constant/options/chfix")) { 
cellSet selectedChfix(mesh, "chfix");cellsChfix0 = selectedChfix.toc();
while (time>=tChfix) {output = getCbuffer("chfix",itChfix,cellsChfix0.size());itChfix +=1;tChfix = output.t;cellsChfixData=output.d;}
cellsChfix = indexC(cellsChfix0,cellsChfixData) ;
}
if (fexists(cur_dir+"/constant/options/cfix")) { 
cellSet selectedCfix(mesh, "cfix");cellsCfix0 = selectedCfix.toc();
while (time>=tCfix) {output = getCbuffer("cfix",itCfix,cellsCfix0.size());itCfix +=1;tCfix = output.t;cellsCfixData=output.d;}
cellsCfix = indexC(cellsCfix0,cellsCfixData) ;
}
if (fexists(cur_dir+"/constant/options/cwel")) { 
cellSet selectedCwel(mesh, "cwel");cellsCwel0 = selectedCwel.toc();
while (time>=tCwel) {output = getCbuffer("cwel",itCwel,cellsCwel0.size());itCwel +=1;tCwel = output.t;cellsCwelData=output.d;}
cellsCwel = indexC(cellsCwel0,cellsCwelData) ;
}
if (fexists(cur_dir+"/constant/options/cghb")) { 
cellSet selectedCghb(mesh, "cghb");cellsCghb0 = selectedCghb.toc();
while (time>=tCghb) {output = getCbuffer("cghb",itCghb,cellsCghb0.size());itCghb +=1;tCghb = output.t;cellsCghbData=output.d;}
cellsCghb = indexC(cellsCghb0,cellsCghbData) ;
}
if (fexists(cur_dir+"/constant/options/criv")) { 
cellSet selectedCriv(mesh, "criv");cellsCriv0 = selectedCriv.toc();
while (time>=tCriv) {output = getCbuffer("criv",itCriv,cellsCriv0.size());itCriv +=1;tCriv = output.t;cellsCrivData=output.d;}
cellsCriv = indexC(cellsCriv0,cellsCrivData) ;
}
if (fexists(cur_dir+"/constant/options/crch")) { 
cellSet selectedCrch(mesh, "crch");cellsCrch0 = selectedCrch.toc();
while (time>=tCrch) {output = getCbuffer("crch",itCrch,cellsCrch0.size());itCrch +=1;tCrch = output.t;cellsCrchData=output.d;}
cellsCrch = indexC(cellsCrch0,cellsCrchData) ;
}
//get the data



// mass balance
if (activateTransport==1) {std::ofstream fich1(cur_dir+"/observation/budgetC.txt");
fich1 << "time InCfix OutCfix InCwel OutCwel InCdrn OutCdrn InCghb OutCghb InCrch OutCrch ";fich1.close();}

float sqICrch,sqICghb,sqICdrn,sqICfix,sqICwel;
float sqOCrch,sqOCghb,sqOCdrn,sqOCfix,sqOCwel;


