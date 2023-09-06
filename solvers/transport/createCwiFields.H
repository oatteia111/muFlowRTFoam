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
labelList cellsSfx0; std::vector<float> cellsSfxData;std::vector<float> condSfx;
labelList cellsSwl0; std::vector<float> cellsSwlData;std::vector<float> condSwl;
labelList cellsSghb0; std::vector<float> cellsSghbData;std::vector<float> condSghb;
labelList cellsSrch0; std::vector<float> cellsSrchData;std::vector<float> condSrch;
scalar imin=0;

// reading the data for fix
std::ifstream inputSfxdata{cur_dir+"/constant/options/sfix", std::ios::binary}; //
//cellsSfxData = {std::istream_iterator<float>{inputSfxdata}, {}};
while (inputSfxdata.read(reinterpret_cast<char*>(&f0), sizeof(float))) {cellsSfxData.push_back(f0);}
if (cellsSfxData.size()>0) { cellSet selectedSfx(mesh, "sfix"); cellsSfx0 = selectedSfx.toc();}
condSfx.resize(cellsSfx0.size());
std::vector<int> cellsSfx = indexC(cellsSfx0,cellsSfxData) ;

// reading the well data for wells
std::ifstream inputSwldata{cur_dir+"/constant/options/swel", std::ios::binary}; //
//cellsSwlData = {std::istream_iterator<float>{inputSwldata}, {}};
while (inputSwldata.read(reinterpret_cast<char*>(&f0), sizeof(float))) {cellsSwlData.push_back(f0);}
if (cellsSwlData.size()>0) {	cellSet selectedSwl(mesh, "swel");cellsSwl0 = selectedSwl.toc();}
condSwl.resize(cellsSwl0.size());
std::vector<int> cellsSwl = indexC(cellsSwl0,cellsSwlData) ;

// reading the data for ghb
std::ifstream inputSghbdata{cur_dir+"/constant/options/sghb", std::ios::binary}; //
//cellsSghbData = {std::istream_iterator<float>{inputSghbdata}, {}};
while (inputSghbdata.read(reinterpret_cast<char*>(&f0), sizeof(float))) {cellsSghbData.push_back(f0);}
if (cellsSghbData.size()>0) {	cellSet selectedSghb(mesh, "sghb");cellsSghb0 = selectedSghb.toc();}
condSghb.resize(cellsSghb0.size());
std::vector<int> cellsSghb = indexC(cellsSghb0,cellsSghbData) ;

// reading for recharge
std::ifstream inputSrchdata{cur_dir+"/constant/options/srch", std::ios::binary}; //
//cellsSrchData = {std::istream_iterator<float>{inputSrchdata}, {}};
while (inputSrchdata.read(reinterpret_cast<char*>(&f0), sizeof(float))) {cellsSrchData.push_back(f0);}
if (cellsSrchData.size()>0) {	cellSet selectedSrch(mesh, "srch");cellsSrch0 = selectedSrch.toc();}
condSrch.resize(cellsSrch0.size());
clayRch.resize(cellsSrch0.size());
std::vector<int> cellsSrch = indexC(cellsSrch0,cellsSrchData) ;

Info <<"end Cwifields " << endl;