dimensionedScalar cpw(transportProperties.lookupOrDefault("cpw",dimensionedScalar("cpw",dimensionSet(0,2,-2,-1,0,0,0),4182)));
dimensionedScalar lbdaTw(transportProperties.lookupOrDefault("lbdaTw",dimensionedScalar("lbdaTw",dimensionSet(1,1,-3,-1,0,0,0),0.6)));
scalar activateThermal(transportProperties.lookupOrDefault<scalar>("activateThermal",0));
scalar dT=1e-9;scalar dT1=1e-9;

//NB (dimensionType MASS , LENGTH , TIME , TEMPERATURE ,MOLES , CURRENT , LUMINOUS_INTENSITY)
Info<< "Reading field T\n" << endl;
volScalarField T  // Temperature
(
    IOobject("T",runTime.timeName(),mesh,IOobject::READ_IF_PRESENT,IOobject::AUTO_WRITE),mesh,dimensionedScalar("T",dimTemperature,25.)
);


volScalarField cps  // heat capacity of the solid [L2 T-2 °C-1]
(
	IOobject("cps",runTime.constant(),mesh,IOobject::READ_IF_PRESENT,IOobject::AUTO_WRITE),mesh,dimensionedScalar("cps",dimensionSet(0,2,-2,-1,0,0,0),840)
);

volScalarField lbdaTs  // Thermal conductivity of the solid [M L T-3 °C-1]
(
	IOobject("lbdaTs",runTime.constant(),mesh,IOobject::READ_IF_PRESENT,IOobject::AUTO_WRITE),mesh,dimensionedScalar("lbdaTs",dimensionSet(1,1,-3,-1,0,0,0),3.5)
);

volScalarField lbdaT = sw*eps*lbdaTw + (1-eps)*lbdaTs;


// creating the internal Su and Sp fields to input temperature conditions inside the domain
volScalarField::Internal TSu
(
   IOobject("TSu",runTime.timeName(),mesh,IOobject::NO_READ,IOobject::NO_WRITE),
		mesh,
		dimensionedScalar("zero",dimensionSet(1,-3,-1,1,0,0,0),0.0) // the overall eq in in K.s-1 but rhow is included in Tsu
		//dimensioned<scalar>("zero",eqn.dimensions()/dimVolume/psi.dimensions(),0.0), false
);

volScalarField::Internal TSp
(
   IOobject("TSp",runTime.timeName(),mesh,IOobject::NO_READ,IOobject::NO_WRITE),
		mesh,
		dimensionedScalar("zero",dimensionSet(1,-3,-1,0,0,0,0),0.0) 
		//dimensioned<scalar>("zero",eqn.dimensions()/dimVolume/psi.dimensions(),0.0), false
);

//////////////////////// read conditions //////////////////////////
labelList cellsThfix0; std::vector<float> condThfix;std::vector<float> cellsThfixData;float tThfix=0;
labelList cellsTfix0; std::vector<float> condTfix;std::vector<float> cellsTfixData;float tTfix=0;
labelList cellsTwel0; std::vector<float> condTwel;std::vector<float> cellsTwelData;float tTwel=0;
labelList cellsTghb0; std::vector<float> condTghb;std::vector<float> cellsTghbData;float tTghb=0;
labelList cellsTriv0; std::vector<float> condTriv;std::vector<float> cellsTrivData;float tTriv=0;
labelList cellsTrch0; std::vector<float> condTrch;std::vector<float> cellsTrchData;float tTrch=0;
int itThfix,itTfix,itTwel,itTghb,itTrch,itTriv;

std::vector<int> cellsThfix,cellsTfix,cellsTwel,cellsTghb,cellsTriv,cellsTrch;

if (fexists(cur_dir+"/constant/options/thfix")) { 
cellSet selectedThfix(mesh, "thfix");cellsThfix0 = selectedThfix.toc();
while (time>=tThfix) 
{output = getCbuffer("thfix",itThfix,cellsThfix0.size());itThfix +=1;tThfix = output.t;cellsThfixData=output.d;}
cellsThfix = indexC(cellsThfix0,cellsThfixData) ;
}
if (fexists(cur_dir+"/constant/options/tfix")) { 
cellSet selectedTfix(mesh, "tfix");cellsTfix0 = selectedTfix.toc();
while (time>=tTfix) 
{output = getCbuffer("tfix",itTfix,cellsTfix0.size());itTfix +=1;tTfix = output.t;cellsTfixData=output.d;}
cellsTfix = indexC(cellsTfix0,cellsTfixData) ;
}
if (fexists(cur_dir+"/constant/options/twel")) { 
cellSet selectedTwel(mesh, "twel");cellsTwel0 = selectedTwel.toc();
while (time>=tTwel) 
{output = getCbuffer("twel",itTwel,cellsTwel0.size());itTwel +=1;tTwel = output.t;cellsTwelData=output.d;}
cellsTwel = indexC(cellsTwel0,cellsTwelData) ;
}
if (fexists(cur_dir+"/constant/options/tghb")) { 
cellSet selectedTghb(mesh, "tghb");cellsTghb0 = selectedTghb.toc();
while (time>=tTghb) 
{output = getCbuffer("tghb",itTghb,cellsTghb0.size());itTghb +=1;tTghb = output.t;cellsTghbData=output.d;}
cellsTghb = indexC(cellsTghb0,cellsTghbData) ;
}
if (fexists(cur_dir+"/constant/options/triv")) { 
cellSet selectedTriv(mesh, "triv");cellsTriv0 = selectedTriv.toc();
while (time>=tTriv) 
{output = getCbuffer("triv",itTriv,cellsTriv0.size());itTriv +=1;tTriv = output.t;cellsTrivData=output.d;}
cellsTriv = indexC(cellsTriv0,cellsTrivData) ;
}
if (fexists(cur_dir+"/constant/options/trch")) { 
cellSet selectedTrch(mesh, "trch");cellsTrch0 = selectedTrch.toc();
while (time>=tTrch) 
{output = getCbuffer("trch",itTrch,cellsTrch0.size());itTrch +=1;tTrch = output.t;cellsTrchData=output.d;}
cellsTrch = indexC(cellsTrch0,cellsTrchData) ;
}

// mass balance
if (activateThermal==1) {std::ofstream fich1(cur_dir+"/observation/budgetT.txt");
fich1 << "time InTfix OutTfix InTwel OutTwel InTdrn OutTdrn InTghb OutTghb InTrch OutTrch ";
fich1.close();}
float sqITrch,sqITghb,sqITdrn,sqITfix,sqITwel;
float sqOTrch,sqOTghb,sqOTdrn,sqOTfix,sqOTwel;

Info<< "All T read\n" << endl;
