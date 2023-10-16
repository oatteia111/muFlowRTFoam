//#include "inletOutletFvPatchScalarField.H" does not exist
//#include "inletOutletFvPatchFields.H"

/*
dimensionedScalar alphaL(transportProperties.lookupOrDefault("alphaL",dimensionedScalar("alphaL",dimLength,0.0)));
dimensionedScalar alphaT(transportProperties.lookupOrDefault("alphaT",dimensionedScalar("alphaT",dimLength,0.0)));
dimensionedScalar lbdaw(transportProperties.lookupOrDefault("lbdaw",dimensionedScalar("lbdaw",dimless/dimTime,0.0)));
dimensionedScalar Dw0(transportProperties.lookupOrDefault("Dw0",dimensionedScalar("Dw0",dimArea/dimTime,1e-10)));

Info<< "Reading field Deffw if present\n" << endl;
volTensorField Deffw
(
    IOobject
    ("Deffw",runTime.timeName(),mesh,IOobject::READ_IF_PRESENT,IOobject::NO_WRITE),
    mesh,
	dimensionedTensor("",dimArea/dimTime,Zero) 
);

// below from https://www.cfd-online.com/Forums/openfoam-programming-development/113261-multiple-volscalarfields-using-list-volscalarfield.html
Info <<"in Cifield nxyz "<<nxyz<< " n comp "<<ph_ncomp << endl;

// Initialize BCs list for Cc to zero-gradient
wordList bcTypes
	(
	Kh.boundaryField().size(),
	zeroGradientFvPatchScalarField::typeName // works
	);
	

List<word> bcNames;bcNames.append("bc3");
for (int i=0;i<bcNames.size();i++)
	{
	const label inletID = mesh.boundaryMesh().findPatchID(bcNames[i]);
	bcTypes[inletID] = fixedValueFvPatchScalarField::typeName;
	}
*/

List <autoPtr <volScalarField> > Cw(ph_ncomp);
forAll (Cw, i)
{
		word nameCi ("Cw" + name(i));
		Cw[i].reset
		(
			new volScalarField
			(
			 IOobject(nameCi,runTime.timeName(),mesh,IOobject::READ_IF_PRESENT,IOobject::AUTO_WRITE),
			 mesh,
			 dimensionedScalar("",dimMass/dimVol,0.), //freak.solu_conc[i]
			 bcTypes
			)
		);
}

// creating the internal Su and Sp fields to input temperature conditions inside the domain
volScalarField::Internal CiSu
(
   IOobject("CiSu",runTime.timeName(),mesh,IOobject::NO_READ,IOobject::NO_WRITE),
		mesh,
		dimensionedScalar("zero",dimensionSet(1,-3,-1,0,0,0,0),0.0) 
);

volScalarField::Internal CiSp
(
   IOobject("CiSp",runTime.timeName(),mesh,IOobject::NO_READ,IOobject::NO_WRITE),
		mesh,
		dimensionedScalar("zero",dimensionSet(0,0,-1,0,0,0,0),0.0) 
);

// -------------- conditions for solutions -------------
labelList cellsSfix0; std::vector<float> condSfix;std::vector<float> cellsSfixData;int itSfix=0;float tSfix=0;
labelList cellsSwel0; std::vector<float> condSwel;std::vector<float> cellsSwelData;int itSwel=0;float tSwel=0;
labelList cellsSghb0; std::vector<float> condSghb;std::vector<float> cellsSghbData;int itSghb=0;float tSghb=0;
labelList cellsSrch0; std::vector<float> condSrch;std::vector<float> cellsSrchData;int itSrch=0;float tSrch=0;

if (fexists(cur_dir+"/constant/options/sfix")) {cellSet selectedSfix(mesh, "sfix");cellsSfix0 = selectedSfix.toc();}
if (fexists(cur_dir+"/constant/options/swel")) {cellSet selectedSwel(mesh, "swel");cellsSwel0 = selectedSwel.toc();}
if (fexists(cur_dir+"/constant/options/sghb")) {cellSet selectedSghb(mesh, "sghb");cellsSghb0 = selectedSghb.toc();}
if (fexists(cur_dir+"/constant/options/srch")) {cellSet selectedSrch(mesh, "srch");cellsSrch0 = selectedSrch.toc();}

if (time>=tSfix) {output = getCbuffer("sfix",itSfix,cellsSfix0.size());itSfix+=1;tSfix = output.t;cellsSfixData=output.d;}
if (time>=tSwel) {output = getCbuffer("swel",itSwel,cellsSwel0.size());itSwel+=1;tSwel = output.t;cellsSwelData=output.d;}
if (time>=tSghb) {output = getCbuffer("sghb",itSghb,cellsSghb0.size());itSghb+=1;tSghb = output.t;cellsSghbData=output.d;}
if (time>=tSrch) {output = getCbuffer("srch",itSrch,cellsSrch0.size());itSrch+=1;tSrch = output.t;cellsSrchData=output.d;}

condSfix.resize(cellsSfix0.size());
condSwel.resize(cellsSwel0.size());
condSghb.resize(cellsSghb0.size());
condSrch.resize(cellsSrch0.size());
std::vector<int> cellsSfix = indexC(cellsSfix0,cellsSfixData) ;
std::vector<int> cellsSwel = indexC(cellsSwel0,cellsSwelData) ;
std::vector<int> cellsSghb = indexC(cellsSghb0,cellsSghbData) ;
std::vector<int> cellsSrch = indexC(cellsSrch0,cellsSrchData) ;

Info <<"end Cwifields " << endl;