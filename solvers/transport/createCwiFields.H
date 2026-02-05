//#include "inletOutletFvPatchScalarField.H" does not exist
//#include "inletOutletFvPatchFields.H"

/*
//stuff below is kept jus tin case
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
std::vector<double> mxCv;
int ncc;
if (activateEK) {ncc=freak.nspc;}
else {ncc=ph_ncomp;}

List <autoPtr <volScalarField> > Cw(ncc);
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

// creating the internal Su and Sp fields to input conditions inside the domain
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

DynamicList<label> inactiveFaces;
// -------------- conditions for solutions -------------
labelList cellsShfix0; std::vector<float> condShfix;std::vector<float> cellsShfixData;float tShfix=0;int itShfix=0;
labelList cellsSpfix0; std::vector<float> condSpfix;std::vector<float> cellsSpfixData;float tSpfix=0;int itSpfix=0;
labelList cellsSfix0; std::vector<float> condSfix;std::vector<float> cellsSfixData;float tSfix=0;int itSfix=0;
labelList cellsSwel0; std::vector<float> condSwel;std::vector<float> cellsSwelData;float tSwel=0;int itSwel=0;
labelList cellsSghb0; std::vector<float> condSghb;std::vector<float> cellsSghbData;float tSghb=0;int itSghb=0;
labelList cellsSrch0; std::vector<float> condSrch;std::vector<float> cellsSrchData;float tSrch=0;int itSrch=0;
std::vector<int> bcCwi;

std::string s="/constant/options";
std::vector<int> cellsShfix,cellsSpfix,cellsSfix,cellsSwel,cellsSghb,cellsSrch;
if (fexists(cur_dir+s+"/shfix")) {
cellSet selectedShfix(mesh, "shfix");cellsShfix0 = selectedShfix.toc();condShfix.resize(cellsShfix0.size());
while (time>=tShfix) {output = getCbuffer("shfix",itShfix,cellsShfix0.size());itShfix+=1;tShfix = output.t;cellsShfixData=output.d;}
cellsShfix = indexC(cellsShfix0,cellsShfixData) ;
}
if (fexists(cur_dir+s+"/spfix")) {
cellSet selectedSpfix(mesh, "spfix");cellsSpfix0 = selectedSpfix.toc();condSpfix.resize(cellsSpfix0.size());
while (time>=tSpfix) {output = getCbuffer("spfix",itSpfix,cellsSpfix0.size());itSpfix+=1;tSpfix = output.t;cellsSpfixData=output.d;}
cellsSpfix = indexC(cellsSpfix0,cellsSpfixData) ;
}
if (fexists(cur_dir+s+"/sfix")) {
cellSet selectedSfix(mesh, "sfix");cellsSfix0 = selectedSfix.toc();condSfix.resize(cellsSfix0.size());
while (time>=tSfix) {output = getCbuffer("sfix",itSfix,cellsSfix0.size());itSfix+=1;tSfix = output.t;cellsSfixData=output.d;}
cellsSfix = indexC(cellsSfix0,cellsSfixData) ;
}
if (fexists(cur_dir+s+"/swel")) {
cellSet selectedSwel(mesh, "swel");cellsSwel0 = selectedSwel.toc();condSwel.resize(cellsSwel0.size());
while (time>=tSwel) {output = getCbuffer("swel",itSwel,cellsSwel0.size());itSwel+=1;tSwel = output.t;cellsSwelData=output.d;}
cellsSwel = indexC(cellsSwel0,cellsSwelData) ;
}
if (fexists(cur_dir+s+"/sghb")) {
cellSet selectedSghb(mesh, "sghb");cellsSghb0 = selectedSghb.toc();condSghb.resize(cellsSghb0.size());
while (time>=tSghb) {output = getCbuffer("sghb",itSghb,cellsSghb0.size());itSghb+=1;tSghb = output.t;cellsSghbData=output.d;}
cellsSghb = indexC(cellsSghb0,cellsSghbData) ;
}
if (fexists(cur_dir+s+"/srch")) {
cellSet selectedSrch(mesh, "srch");cellsSrch0 = selectedSrch.toc();condSrch.resize(cellsSrch0.size());
while (time>=tSrch) {output = getCbuffer("srch",itSrch,cellsSrch0.size());itSrch+=1;tSrch = output.t;cellsSrchData=output.d;}
cellsSrch = indexC(cellsSrch0,cellsSrchData) ;
}

// for mass balance
std::vector<float> sqISfix(ph_ncomp),sqISwel(ph_ncomp),sqISghb(ph_ncomp),sqISdrn(ph_ncomp),sqISrch(ph_ncomp);
std::vector<float> sqOSfix(ph_ncomp),sqOSwel(ph_ncomp),sqOSghb(ph_ncomp),sqOSdrn(ph_ncomp),sqOSrch(ph_ncomp);
	for (int ic=4;ic<ph_ncomp;ic++)
	{
	std::string s=std::to_string(ic);
	std::ofstream fich1(cur_dir+"/observation/budgetCw"+s+".txt");
	fich1 << "time InCfix OutCfix InCwel OutCwel InCdrn OutCdrn InCghb OutCghb InCrch OutCrch \n";
	fich1.close();
	}


Info <<"end Cwifields \n" << endl;