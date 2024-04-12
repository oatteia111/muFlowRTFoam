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
labelList cellsShfix0; std::vector<float> condShfix;std::vector<float> cellsShfixData;float tShfix=0;int itShfix=0;
labelList cellsSfix0; std::vector<float> condSfix;std::vector<float> cellsSfixData;float tSfix=0;int itSfix=0;
labelList cellsSwel0; std::vector<float> condSwel;std::vector<float> cellsSwelData;float tSwel=0;int itSwel=0;
labelList cellsSghb0; std::vector<float> condSghb;std::vector<float> cellsSghbData;float tSghb=0;int itSghb=0;
labelList cellsSrch0; std::vector<float> condSrch;std::vector<float> cellsSrchData;float tSrch=0;int itSrch=0;

std::string s="/constant/options";
if (fexists(cur_dir+s+"/shfix")) {cellSet selectedShfix(mesh, "shfix");cellsShfix0 = selectedShfix.toc();}
if (fexists(cur_dir+s+"/sfix")) {cellSet selectedSfix(mesh, "sfix");cellsSfix0 = selectedSfix.toc();}
if (fexists(cur_dir+s+"/swel")) {cellSet selectedSwel(mesh, "swel");cellsSwel0 = selectedSwel.toc();}
if (fexists(cur_dir+s+"/sghb")) {cellSet selectedSghb(mesh, "sghb");cellsSghb0 = selectedSghb.toc();}
if (fexists(cur_dir+s+"/srch")) {cellSet selectedSrch(mesh, "srch");cellsSrch0 = selectedSrch.toc();}

time = mesh.time().value()/day;
while (time>=tShfix) {output = getCbuffer("shfix",itShfix,cellsShfix0.size());itShfix+=1;tShfix = output.t;cellsShfixData=output.d;tnext=tShfix;}
while (time>=tSfix) {output = getCbuffer("sfix",itSfix,cellsSfix0.size());itSfix+=1;tSfix = output.t;cellsSfixData=output.d;tnext=tSfix;}
while (time>=tSwel) {output = getCbuffer("swel",itSwel,cellsSwel0.size());itSwel+=1;tSwel = output.t;cellsSwelData=output.d;tnext=tSwel;}
while (time>=tSghb) {output = getCbuffer("sghb",itSghb,cellsSghb0.size());itSghb+=1;tSghb = output.t;cellsSghbData=output.d;tnext=tSghb;}
while (time>=tSrch) {output = getCbuffer("srch",itSrch,cellsSrch0.size());itSrch+=1;tSrch = output.t;cellsSrchData=output.d;tnext=tSrch;}

condShfix.resize(cellsShfix0.size());
condSfix.resize(cellsSfix0.size());
condSwel.resize(cellsSwel0.size());
condSghb.resize(cellsSghb0.size());
condSrch.resize(cellsSrch0.size());
std::vector<int> cellsShfix = indexC(cellsShfix0,cellsShfixData) ;
std::vector<int> cellsSfix = indexC(cellsSfix0,cellsSfixData) ;
std::vector<int> cellsSwel = indexC(cellsSwel0,cellsSwelData) ;
std::vector<int> cellsSghb = indexC(cellsSghb0,cellsSghbData) ;
std::vector<int> cellsSrch = indexC(cellsSrch0,cellsSrchData) ;

// for mass balance
std::vector<float> sqISfix(ph_ncomp),sqISwel(ph_ncomp),sqISghb(ph_ncomp),sqISdrn(ph_ncomp),sqISrch(ph_ncomp);
std::vector<float> sqOSfix(ph_ncomp),sqOSwel(ph_ncomp),sqOSghb(ph_ncomp),sqOSdrn(ph_ncomp),sqOSrch(ph_ncomp);
	for (int ic=0;ic<ph_ncomp;ic++)
	{
	std::string s=std::to_string(ic);
	std::ofstream fich1(cur_dir+"/observation/budgetCw"+s+".txt");// Ouvrir le fichier en mode append
	fich1 <<"t sIFix sOFix sIWel sOWel sIDrn sODrn sIGhb sOGhb"<< std::endl;
	fich1.close();
	}


Info <<"end Cwifields " << endl;