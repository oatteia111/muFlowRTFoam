dimensionedScalar Dg0(transportProperties.lookupOrDefault("Dg0",dimensionedScalar("Dg0",dimArea/dimTime,1e-6)));
scalar diffusionEqn(transportProperties.lookupOrDefault<scalar>("diffusionEqn",0));

Info<< "Reading field Deffg if present\n" << endl;
volTensorField Deffg
(
    IOobject
    ("Deffg",runTime.timeName(),mesh,IOobject::READ_IF_PRESENT,IOobject::NO_WRITE),
    mesh,
	dimensionedTensor("",dimArea/dimTime,Zero) 
);

Info <<"in Cgfield nxyz "<<nxyz<< " g comp "<<ph_gcomp << endl;

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

// -------------- conditions for gas mixtures -------------
labelList cellsGfix0; std::vector<float> qGfix;std::vector<float> cellsGfixData;std::vector<float> condGfix;
labelList cellsGwel0; std::vector<float> qGwel;std::vector<float> cellsGwelData;std::vector<float> condGwel;
int itGfix=0;int itGwel=0;float tGfix=0;float tGwel=0;

if (fexists(cur_dir+"/constant/options/gfix")) { cellSet selectedGfix(mesh, "gfix");cellsGfix0 = selectedGfix.toc();}
if (fexists(cur_dir+"/constant/options/gwel")) { cellSet selectedGwel(mesh, "gwel");cellsGwel0 = selectedGwel.toc();}

while (time>=tGfix) {output = getCbuffer("gfix",itGfix,cellsGfix0.size());itGfix+=1;tGfix = output.t;cellsGfixData=output.d;}
while (time>=tGwel) {output = getCbuffer("gwel",itGwel,cellsGwel0.size());itGwel+=1;tGwel = output.t;cellsGwelData=output.d;}

qGfix.resize(cellsGfix0.size());
qGwel.resize(cellsGwel0.size());
condGfix.resize(cellsGfix0.size());
condGwel.resize(cellsGwel0.size());

std::vector<int> cellsGfix = indexC(cellsGfix0,cellsGfixData) ;
std::vector<int> cellsGwel = indexC(cellsGwel0,cellsGwelData) ;

Info <<"end Cgifields " << endl;