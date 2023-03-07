// transprot properties variable
Info << "Reading transportProperties\n" << endl;
IOdictionary transportProperties
(
    IOobject
    ("transportProperties",runTime.constant(), mesh,IOobject::MUST_READ_IF_MODIFIED,IOobject::NO_WRITE)
);
scalar flowStartSteady(transportProperties.lookupOrDefault<scalar>("flowStartSteady",0.));
scalar flowType(transportProperties.lookupOrDefault<scalar>("flowType",1)); //0 no flow, 1 confined, 2 unconfined, 3 unsat 4 dualphase
scalar activateCapillarity(transportProperties.lookupOrDefault<scalar>("activateCapillarity",1));
scalar activateTransport(transportProperties.lookupOrDefault<scalar>("activateTransport",0));
scalar activateReaction(transportProperties.lookupOrDefault<scalar>("activateReaction",0));
scalar reactionSteps(transportProperties.lookupOrDefault<scalar>("reactionSteps",10));
//scalar phreeqcVm(transportProperties.lookupOrDefault<scalar>("phreeqcVm",24.5));
//scalar phreeqcP(transportProperties.lookupOrDefault<scalar>("phreeqcP",1.));
scalar stor0(transportProperties.lookupOrDefault<scalar>("stor0",1e-4));
int nlay(transportProperties.lookupOrDefault<int>("nlay",1));
int ncell_lay(transportProperties.lookupOrDefault<int>("ncell_lay",1));
int ncell=nlay*ncell_lay;
Info << "transportP read "<< endl;

Info << nl << "Reading field h\n" << endl;
volScalarField eps
(
    IOobject("eps",runTime.constant(),mesh,IOobject::READ_IF_PRESENT,IOobject::NO_WRITE),mesh,dimensionedScalar("eps",dimless,0.3)
);
volScalarField h
(
    IOobject("h",runTime.timeName(),mesh,IOobject::READ_IF_PRESENT,IOobject::AUTO_WRITE),mesh,dimensionedScalar("h",dimLength,0)
);
volScalarField hp
(
    IOobject("hp",runTime.timeName(),mesh,IOobject::READ_IF_PRESENT,IOobject::AUTO_WRITE),mesh,dimensionedScalar("hp",dimLength,0)
);
volScalarField zbot
(
    IOobject("zbot",runTime.constant(),mesh,IOobject::READ_IF_PRESENT,IOobject::NO_WRITE),mesh,dimensionedScalar("zbot",dimLength,0.)
);
volScalarField thk
(
    IOobject("thk",runTime.constant(),mesh,IOobject::READ_IF_PRESENT,IOobject::NO_WRITE),mesh,dimensionedScalar("thk",dimLength,1)
);
volScalarField p
(
    IOobject("p",runTime.timeName(),mesh,IOobject::READ_IF_PRESENT,IOobject::AUTO_WRITE),mesh,dimensionedScalar("p",dimensionSet(1,-1,-2,0,0,0,0),0)
);
volScalarField sw
(
    IOobject("sw",runTime.timeName(),mesh,IOobject::READ_IF_PRESENT,IOobject::AUTO_WRITE),mesh,dimensionedScalar("sw",dimless,1.)
);
volScalarField phreeqcVm
(
    IOobject("phreeqcVm",runTime.constant(),mesh,IOobject::READ_IF_PRESENT,IOobject::NO_WRITE),mesh,dimensionedScalar("phreeqcVm",dimless,24.5)
);

volScalarField unity
(
    IOobject("unity",runTime.constant(),mesh,IOobject::NO_READ,IOobject::NO_WRITE),mesh,dimensionedScalar("unity",dimless,1)
);
volScalarField a("a",unity);
a.dimensions().reset(dimMass/dimVol); //for C
volScalarField b("b",unity);
b.dimensions().reset(dimless/dimLength); //for fvOptions

//Info<<"eps "<<eps[0]<<endl;
// time step data
const bool adjustTimeStep = runTime.controlDict().lookupOrDefault("adjustTimeStep", true);
const word timeStepControl = runTime.controlDict().lookupOrDefault<word>("timeStepControl", "");
word CFL =  runTime.controlDict().lookupOrDefault<word>("CFL", "Coats");
//if (adjustTimeStep) {adjustTimeStep=true;} // to remove warnings at compilation
int CFLint = -1;
if (CFL == "Coats") CFLint = 0; 
else if (CFL == "Courant") CFLint = 1; 
else{ FatalErrorIn ("in createFields.H") << "CFL condition unrecongnized : Coats and Courant available" << exit(FatalError);}
Info<<"CFL "<<CFL<<" int "<<CFLint<<endl;

volScalarField stor("stor", unity); //Info << "stor " << stor().size() << endl;
forAll(stor(), celli) { stor[celli] = stor0/thk[celli] ;}

//phases 
autoPtr<incompressiblePhase> phasew = incompressiblePhase::New(mesh,transportProperties,"w");
volVectorField& Uw = phasew->U();
const dimensionedScalar& rhow = phasew->rho();
const dimensionedScalar& muw = phasew->mu();
phasew->phi().writeOpt()=IOobject::NO_WRITE;
Info << "phase done "<< endl;
volVectorField& U = Uw;
surfaceScalarField& phiw = phasew->phi();

volScalarField dffPhi("dffPhi",unity); //flow balance 
dffPhi.dimensions().reset(dimVol/dimTime); 

/*surfaceScalarField phiw
(
    IOobject("phiw", runTime.timeName(), mesh,IOobject::READ_IF_PRESENT,IOobject::AUTO_WRITE),mesh,dimensionedScalar("phiw",dimVol/dimTime,0.)
);*/

Info<< "reading Kh and Kv field " << endl;

volScalarField Kh
(
    IOobject("Kh",runTime.constant(),mesh,IOobject::MUST_READ,IOobject::AUTO_WRITE),mesh
);
volScalarField Kv
(
    IOobject("Kv",runTime.constant(),mesh,IOobject::MUST_READ,IOobject::AUTO_WRITE),mesh
);

surfaceScalarField Kf
(
    IOobject("Kf",runTime.timeName(),mesh,IOobject::NO_READ,IOobject::NO_WRITE),mesh,dimensionedScalar("Kf",dimArea,1.e-12)
);
Kf = fvc::interpolate(Kh,"Kh");

//create a surface scalar field with harmonic mean different for the vertical and non-vertical faces
label oCell , nCell;
scalar odCface,ndCface;
const surfaceVectorField& oCf  = mesh.Cf();
const volVectorField& oC = mesh.C();
forAll( mesh.owner() , iface ) 
{
	// Get the left and right cell index
	oCell = mesh.owner()[iface];odCface = mag(oCf[iface]-oC[oCell]);
	nCell = mesh.neighbour()[iface];ndCface = mag(oCf[iface]-oC[nCell]);
	if (abs(mesh.Sf()[iface].z())<mag(mesh.Sf()[iface])*1e-5)  //vertical faces
	{
	//Kf[iface] = 2*(Kh[oCell] * Kh[nCell])/(Kh[oCell] + Kh[nCell]); // harmonic mean in x direction
	Kf[iface] = (Kh[oCell]*Kh[nCell]) /(odCface*Kh[oCell] + ndCface*Kh[nCell] )*(ndCface+odCface);//mag(mesh.Sf()[iface])*
	if (iface==10) {Info<<"Kf 10 "<<Kf[10]<<endl;}
	}
	else
	{
	//Kf[iface] = 2*(Kv[oCell] * Kv[nCell])/(Kv[oCell] + Kv[nCell]); // harmonic mean in z direction
	Kf[iface] = (Kv[oCell]*Kv[nCell]) /(odCface*Kv[oCell] + ndCface*Kv[nCell] )*(ndCface+odCface); // Lown*Kneighb normal
	}
	//Info << " ifac "<<iface << " K "<<Kf[iface] <<endl;
}

//set values for sw
//forAll(swT, celli) { swT[celli] = tensor::I;}//.xx() = 1; swT[celli].yy() = 1;swT[celli].zz() = 1;}
surfaceScalarField M0f ("M0f",rhow*mag(g)*Kf/muw);
surfaceScalarField Mf ("Mf",M0f*1);
Info<< "end Kf tensor \n" << endl;

//various
int i,j; int icount = 0;float day = 86400.;
std::vector<double> condHw;

// -------------- conditions for head -------------
// reading the hwell data for hfix
labelList cellsHfx;
std::vector<float> cellsHfxData;
cellSet selectedCells(mesh, "hfix");
cellsHfx = selectedCells.toc();
std::ifstream inputHfxdata{cur_dir+"/constant/options/hfix"}; //
cellsHfxData = {std::istream_iterator<float>{inputHfxdata}, {}};

// reading the hwell data for wells
labelList cellsHw;
std::vector<float> cellsHwData;
cellSet selectedCells1(mesh, "hwel");
cellsHw = selectedCells1.toc();
std::ifstream inputHwdata{cur_dir+"/constant/options/hwel"}; //
cellsHwData = {std::istream_iterator<float>{inputHwdata}, {}};

