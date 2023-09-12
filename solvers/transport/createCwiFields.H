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
labelList cellsSfix0; std::vector<float> condSfix;
labelList cellsSwel0; std::vector<float> condSwel;
labelList cellsSghb0; std::vector<float> condSghb;
labelList cellsSrch0; std::vector<float> condSrch;

if (fexists(cur_dir+"/constant/options/sfix")) {cellSet selectedSfix(mesh, "sfix");cellsSfix0 = selectedSfix.toc();}
if (fexists(cur_dir+"/constant/options/swel")) {cellSet selectedSwel(mesh, "swel");cellsSwel0 = selectedSwel.toc();}
if (fexists(cur_dir+"/constant/options/sghb")) {cellSet selectedSghb(mesh, "sghb");cellsSghb0 = selectedSghb.toc();}
if (fexists(cur_dir+"/constant/options/srch")) {cellSet selectedSrch(mesh, "srch");cellsSrch0 = selectedSrch.toc();}

int itSfix=0;output = getCbuffer("sfix",itSfix,cellsSfix0.size()); 
float tSfix = output.t; std::vector<float> cellsSfixData=output.d;

int itSwel=0;output = getCbuffer("swel",itSwel,cellsSwel0.size()); 
float tSwel = output.t; std::vector<float> cellsSwelData=output.d;

int itSghb=0;output = getCbuffer("sghb",itSghb,cellsSghb0.size()); 
float tSghb = output.t; std::vector<float> cellsSghbData=output.d;

int itSrch=0;output = getCbuffer("srch",itSrch,cellsSrch0.size()); 
float tSrch = output.t; std::vector<float> cellsSrchData=output.d;

condSfix.resize(cellsSfix0.size());
condSwel.resize(cellsSwel0.size());
condSghb.resize(cellsSghb0.size());
condSrch.resize(cellsSrch0.size());
std::vector<int> cellsSfix = indexC(cellsSfix0,cellsSfixData) ;
std::vector<int> cellsSwel = indexC(cellsSwel0,cellsSwelData) ;
std::vector<int> cellsSghb = indexC(cellsSghb0,cellsSghbData) ;
std::vector<int> cellsSrch = indexC(cellsSrch0,cellsSrchData) ;
/*
std::vector<float> cellsSfxData;std::vector<float> condSfx;
labelList cellsSwl0; std::vector<float> cellsSwlData;std::vector<float> condSwl;
labelList cellsSghb0; std::vector<float> cellsSghbData;std::vector<float> condSghb;
labelList cellsSrch0; std::vector<float> cellsSrchData;std::vector<float> condSrch;
scalar imin=0;


// reading the well data for wells
std::ifstream inputSwldata{cur_dir+"/constant/options/swel", std::ios::binary}; //
//cellsSwlData = {std::istream_iterator<float>{inputSwldata}, {}};
while (inputSwldata.read(reinterpret_cast<char*>(&f0), sizeof(float))) {cellsSwlData.push_back(f0);}
if (cellsSwlData.size()>0) {	cellSet selectedSwl(mesh, "swel");cellsSwl0 = selectedSwl.toc();}

// reading the data for ghb
std::ifstream inputSghbdata{cur_dir+"/constant/options/sghb", std::ios::binary}; //
//cellsSghbData = {std::istream_iterator<float>{inputSghbdata}, {}};
while (inputSghbdata.read(reinterpret_cast<char*>(&f0), sizeof(float))) {cellsSghbData.push_back(f0);}
if (cellsSghbData.size()>0) {	cellSet selectedSghb(mesh, "sghb");cellsSghb0 = selectedSghb.toc();}


// reading for recharge
std::ifstream inputSrchdata{cur_dir+"/constant/options/srch", std::ios::binary}; //
//cellsSrchData = {std::istream_iterator<float>{inputSrchdata}, {}};
while (inputSrchdata.read(reinterpret_cast<char*>(&f0), sizeof(float))) {cellsSrchData.push_back(f0);}
if (cellsSrchData.size()>0) {	cellSet selectedSrch(mesh, "srch");cellsSrch0 = selectedSrch.toc();}
clayRch.resize(cellsSrch0.size());
*/

Info <<"end Cwifields " << endl;