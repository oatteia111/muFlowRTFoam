dimensionedScalar lbdag(transportProperties.lookupOrDefault("lbdag",dimensionedScalar("lbdag",dimless/dimTime,0.0)));
dimensionedScalar Dg0(transportProperties.lookupOrDefault("Dg0",dimensionedScalar("Dg0",dimArea/dimTime,1e-6)));
scalar activateEbullition(transportProperties.lookupOrDefault("activateEbullition",0));
volScalarField vDg0(IOobject("vDg0",runTime.constant(),mesh,IOobject::READ_IF_PRESENT,IOobject::NO_WRITE),mesh,dimensionedScalar("vDg0",dimArea/dimTime,1.e-12) );
if (vDg0[0]==1e-12) {vDg0 = Dg0;}
scalar diffusionEqn(transportProperties.lookupOrDefault<scalar>("diffusionEqn",0));

Info<< "Reading field factDg if present\n" << endl;
volTensorField Deffg
(IOobject("Deffg",runTime.timeName(),mesh,IOobject::READ_IF_PRESENT,IOobject::AUTO_WRITE),mesh,dimensionedTensor("",dimArea/dimTime,Zero) );

volScalarField factDg( IOobject("factDg",runTime.timeName(),mesh,IOobject::READ_IF_PRESENT,IOobject::NO_WRITE),mesh,dimensionedScalar("",dimless,1.) );

Info <<"in Cgfield nxyz "<<nxyz<< " g comp "<<ph_gcomp << " factDg "<<factDg[0]<<" "<<factDg[1]<<endl;

// Initialize BCs list for Cc to zero-gradient
/*
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
List <autoPtr <volScalarField> > Cg(ph_gcomp);
forAll (Cg, i)
{
		word nameCi ("Cg" + name(i));
		Cg[i].reset
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
volScalarField::Internal GiSu
(
   IOobject("GiSu",runTime.timeName(),mesh,IOobject::NO_READ,IOobject::NO_WRITE),
		mesh,
		dimensionedScalar("zero",dimensionSet(1,-3,-1,0,0,0,0),0.0) 
);

volScalarField::Internal GiSp
(
   IOobject("GiSp",runTime.timeName(),mesh,IOobject::NO_READ,IOobject::NO_WRITE),
		mesh,
		dimensionedScalar("zero",dimensionSet(0,0,-1,0,0,0,0),0.0) 
);

// -------------- conditions for gas mixtures -------------
labelList cellsGfix0; std::vector<float> qGfix;std::vector<float> cellsGfixData;std::vector<float> condGfix;
labelList cellsGwel0; std::vector<float> qGwel;std::vector<float> cellsGwelData;std::vector<float> condGwel;
int itGfix=0;int itGwel=0;float tGfix=0;float tGwel=0;

std::vector<int> cellsGfix,cellsGwel;

if (fexists(cur_dir+"/constant/options/gfix")) { 
cellSet selectedGfix(mesh, "gfix");cellsGfix0 = selectedGfix.toc();qGfix.resize(cellsGfix0.size());condGfix.resize(cellsGfix0.size());
while (time>=tGfix) {output = getCbuffer("gfix",itGfix,cellsGfix0.size());itGfix+=1;tGfix = output.t;cellsGfixData=output.d;}
cellsGfix = indexC(cellsGfix0,cellsGfixData) ;
}
if (fexists(cur_dir+"/constant/options/gwel")) { 
cellSet selectedGwel(mesh, "gwel");cellsGwel0 = selectedGwel.toc();qGwel.resize(cellsGwel0.size());condGwel.resize(cellsGwel0.size());
while (time>=tGwel) {output = getCbuffer("gwel",itGwel,cellsGwel0.size());itGwel+=1;tGwel = output.t;cellsGwelData=output.d;}
cellsGwel = indexC(cellsGwel0,cellsGwelData) ;
}


Info <<"end Cgifields " << endl;