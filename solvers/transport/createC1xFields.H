Info<< "Reading C1x transport props\n" << endl;

scalar Kdg(transportProperties.lookupOrDefault<scalar>("Kdg",0.));
scalar Kds(transportProperties.lookupOrDefault<scalar>("Kds",0.));
dimensionedScalar kmg(transportProperties.lookupOrDefault("kmg",dimensionedScalar("kmg",dimless/dimTime,0.0)));
dimensionedScalar kms(transportProperties.lookupOrDefault("kms",dimensionedScalar("kms",dimless/dimTime,0.0)));
dimensionedScalar lbdas(transportProperties.lookupOrDefault("lbdas",dimensionedScalar("lbdas",dimless/dimTime,0.0)));

volScalarField C1w
(
    IOobject("C1w",runTime.timeName(), mesh,IOobject::READ_IF_PRESENT,IOobject::AUTO_WRITE),mesh,dimensionedScalar("",dimMass/dimVol,0) 
);
volScalarField C1g
(
    IOobject("C1g",runTime.timeName(), mesh,IOobject::READ_IF_PRESENT,IOobject::AUTO_WRITE),mesh,dimensionedScalar("",dimMass/dimVol,0) 
);
volScalarField C1s
(
    IOobject("C1s",runTime.timeName(),mesh,IOobject::READ_IF_PRESENT,IOobject::NO_WRITE),mesh,dimensionedScalar("",dimMass/dimVol,0) 
);

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

// -------------- conditions for concentration -------------
// thre is chfix: conc for the cells where h is fixed and cfix: ccells for which only c is fixed
labelList cellsCwhfix0; std::vector<float> cellsCwhfixData;float tCwhfix=0;int itCwhfix=0;
labelList cellsCwfix0; std::vector<float> cellsCwfixData;float tCwfix=0;int itCwfix=0;
labelList cellsCwwel0; std::vector<float> cellsCwwelData;float tCwwel=0;int itCwwel=0;
labelList cellsCghfix0; std::vector<float> cellsCghfixData;float tCghfix=0;int itCghfix=0;
labelList cellsCgfix0; std::vector<float> cellsCgfixData;float tCgfix=0;int itCgfix=0;
labelList cellsCgwel0; std::vector<float> cellsCgwelData;float tCgwel=0;int itCgwel=0;

// reading the data for hfix
time = mesh.time().value() ;
std::vector<int> cellsCwhfix,cellsCwfix,cellsCwwel,cellsCghfix,cellsCgfix,cellsCgwel;

if (fexists(cur_dir+"/constant/options/cwhfix")) { 
cellSet selectedCwhfix(mesh, "cwhfix");cellsCwhfix0 = selectedCwhfix.toc();
while (time>=tCwhfix) {output = getCbuffer("cwhfix",itChfix,cellsCwhfix0.size());itCwhfix +=1;tCwhfix = output.t;cellsCwhfixData=output.d;}
cellsCwhfix = indexC(cellsCwhfix0,cellsCwhfixData) ;
}
if (fexists(cur_dir+"/constant/options/cwfix")) { 
cellSet selectedCwfix(mesh, "cwfix");cellsCfix0 = selectedCwfix.toc();
while (time>=tCwfix) {output = getCbuffer("cwfix",itCwfix,cellsCwfix0.size());itCwfix +=1;tCwfix = output.t;cellsCwfixData=output.d;}
cellsCwfix = indexC(cellsCwfix0,cellsCwfixData) ;
}
if (fexists(cur_dir+"/constant/options/cwwel")) { 
cellSet selectedCwwel(mesh, "cwwel");cellsCwwel0 = selectedCwwel.toc();
while (time>=tCwwel) {output = getCbuffer("cwwel",itCwwel,cellsCwwel0.size());itCwwel +=1;tCwwel = output.t;cellsCwwelData=output.d;}
cellsCwwel = indexC(cellsCwwel0,cellsCwwelData) ;
}
if (fexists(cur_dir+"/constant/options/cgfix")) { 
cellSet selectedCgfix(mesh, "cgfix");cellsCfix0 = selectedCgfix.toc();
while (time>=tCgfix) {output = getCbuffer("cgfix",itCgfix,cellsCgfix0.size());itCgfix +=1;tCgfix = output.t;cellsCgfixData=output.d;}
cellsCgfix = indexC(cellsCgfix0,cellsCgfixData) ;
}
if (fexists(cur_dir+"/constant/options/cgwel")) { 
cellSet selectedCgwel(mesh, "cgwel");cellsCgwel0 = selectedCgwel.toc();
while (time>=tCgwel) {output = getCbuffer("cgwel",itCgwel,cellsCgwel0.size());itCgwel +=1;tCgwel = output.t;cellsCgwelData=output.d;}
cellsCgwel = indexC(cellsCgwel0,cellsCgwelData) ;
}
