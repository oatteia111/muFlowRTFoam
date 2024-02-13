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
scalar deltaTchem(transportProperties.lookupOrDefault<scalar>("deltaTchem",86400));
scalar stor0(transportProperties.lookupOrDefault<scalar>("stor0",1e-4));
scalar maxHresid(transportProperties.lookupOrDefault<scalar>("maxHresid",1e-4));
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
volScalarField a1("a1",unity);
a1.dimensions().reset(dimless); //for othe ruse
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

/*//phases 
autoPtr<incompressiblePhase> phasew = incompressiblePhase::New(mesh,transportProperties,"w");
volVectorField& Uw = phasew->U();
const dimensionedScalar& rhow = phasew->rho();
const dimensionedScalar& muw = phasew->mu();
phasew->phi().writeOpt()=IOobject::NO_WRITE;
Info << "phase done "<< endl;
volVectorField& U = Uw;
surfaceScalarField& phiw = phasew->phi();*/
dimensionedScalar rhow(transportProperties.lookupOrDefault("rhow",dimensionedScalar("rhow",dimensionSet(1,-3,0,0,0,0,0),1e3)));
dimensionedScalar muw(transportProperties.lookupOrDefault("muw",dimensionedScalar("muw",dimensionSet(1,-1,-1,0,0,0,0),1e-3)));
volVectorField Uw
(
    IOobject("Uw",runTime.timeName(),mesh,IOobject::READ_IF_PRESENT,IOobject::AUTO_WRITE),mesh,dimensionedVector("Uw",dimLength/dimTime,vector::zero)
);
surfaceScalarField phiw
(
    IOobject("phiw", runTime.timeName(), mesh,IOobject::READ_IF_PRESENT,IOobject::AUTO_WRITE),mesh,dimensionedScalar("phiw",dimVol/dimTime,0.)
);

//hydraulic conductivity
Info<< "reading Kh and Kv field " << endl;

volScalarField Kh
(
    IOobject("Kh",runTime.constant(),mesh,IOobject::MUST_READ,IOobject::AUTO_WRITE),mesh
);
volScalarField Kv
(
    IOobject("Kv",runTime.constant(),mesh,IOobject::MUST_READ,IOobject::AUTO_WRITE),mesh
);
volScalarField storC
(
    IOobject("storC",runTime.constant(),mesh,IOobject::READ_IF_PRESENT,IOobject::AUTO_WRITE),mesh,dimensionedScalar("storC",dimless,1.e-6)
);
volScalarField storU
(
    IOobject("storU",runTime.constant(),mesh,IOobject::READ_IF_PRESENT,IOobject::AUTO_WRITE),mesh,dimensionedScalar("storU",dimless,0.25)
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
	//if (iface==10) {Info<<"Kf 10 "<<Kf[10]<<endl;}
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
int i,i1,j,ilay,ilay1; int icount = 0;double day = 86400.;double v=0;
volScalarField carea ("carea",unity); // this is the cell area but unitless
forAll(carea(), celli) { carea[celli] = mesh.V()[celli]/thk[celli] ;}

//deltat
scalar maxDeltaT = runTime.controlDict().lookupOrDefault("maxDeltaT", GREAT);
scalar newDeltaT;

// creating the internal Su and Sp fields to input head conditions inside the domain
volScalarField::Internal hSu
(
   IOobject("hSu",runTime.timeName(),mesh,IOobject::NO_READ,IOobject::NO_WRITE),
		mesh,
		dimensionedScalar("zero",dimensionSet(0,0,-1,0,0,0,0),0.0) // the overall eq in in m3.s-1
);

volScalarField::Internal hSp
(
   IOobject("hSp",runTime.timeName(),mesh,IOobject::NO_READ,IOobject::NO_WRITE),
		mesh,
		dimensionedScalar("zero",dimensionSet(0,-1,-1,0,0,0,0),0.0) 
);


// -------------- conditions for head -------------
labelList cellsHfix0; std::vector<float> qHfix;std::vector<float> cellsHfixData;float tHfix=0;int itHfix=0;
labelList cellsHwel0; std::vector<float> qHwel;std::vector<float> cellsHwelData;float tHwel=0;int itHwel=0;
labelList cellsHdrn0; std::vector<float> qHdrn;std::vector<float> cellsHdrnData;float tHdrn=0;int itHdrn=0;
labelList cellsHriv0; std::vector<float> qHriv;std::vector<float> cellsHrivData;float tHriv=0;int itHriv=0; 
labelList cellsHghb0; std::vector<float> qHghb;std::vector<float> cellsHghbData;float tHghb=0;int itHghb=0; std::vector<float> qGerr;
labelList cellsHrch0; std::vector<float> qHrch;std::vector<float> cellsHrchData;float tHrch=0;int itHrch=0; std::vector<float> clayRch;

// reading the data for hfix
//std::ifstream inputHfxdata{cur_dir+"/constant/options/hfix"}; //
//cellsHfxData = {std::istream_iterator<double>{inputHfxdata}, {}};
if (fexists(cur_dir+"/constant/options/hfix")) { cellSet selectedHfix(mesh, "hfix");cellsHfix0 = selectedHfix.toc();}
if (fexists(cur_dir+"/constant/options/hwel")) { cellSet selectedHwel(mesh, "hwel");cellsHwel0 = selectedHwel.toc();}
if (fexists(cur_dir+"/constant/options/hdrn")) { cellSet selectedHdrn(mesh, "hdrn");cellsHdrn0 = selectedHdrn.toc();}
if (fexists(cur_dir+"/constant/options/hriv")) { cellSet selectedHriv(mesh, "hriv");cellsHriv0 = selectedHriv.toc();}
if (fexists(cur_dir+"/constant/options/hghb")) { cellSet selectedHghb(mesh, "hghb");cellsHghb0 = selectedHghb.toc();}
if (fexists(cur_dir+"/constant/options/hrch")) { cellSet selectedHrch(mesh, "hrch");cellsHrch0 = selectedHrch.toc();}
//get the data
outData output;float time = mesh.time().value()/day;float tnext;//tnext does not work that way, not used
while (time>=tHfix) {output = getCbuffer("hfix",itHfix,cellsHfix0.size());itHfix +=1;tHfix = output.t;cellsHfixData=output.d;}
while (time>=tHwel) {output = getCbuffer("hwel",itHwel,cellsHwel0.size());itHwel +=1;tHwel = output.t;cellsHwelData=output.d;}
while (time>=tHdrn) {output = getCbuffer("hdrn",itHdrn,cellsHdrn0.size());itHdrn +=1;tHdrn = output.t;cellsHdrnData=output.d;}
while (time>=tHriv) {output = getCbuffer("hriv",itHriv,cellsHriv0.size());itHriv +=1;tHriv = output.t;cellsHrivData=output.d;}
while (time>=tHghb) {output = getCbuffer("hghb",itHghb,cellsHghb0.size());itHghb +=1;tHghb = output.t;cellsHghbData=output.d;}
while (time>=tHrch) {output = getCbuffer("hrch",itHrch,cellsHrch0.size());itHrch +=1;tHrch = output.t;cellsHrchData=output.d;}
//tnext should be equal to the min of all t (but not 0)
qHfix.resize(cellsHfix0.size());
qHwel.resize(cellsHwel0.size());
qHdrn.resize(cellsHdrn0.size());
qHriv.resize(cellsHriv0.size());
qHghb.resize(cellsHghb0.size());qGerr.resize(cellsHghb0.size());
qHrch.resize(cellsHrch0.size());clayRch.resize(cellsHrch0.size());
std::vector<int> cellsHfix = indexC(cellsHfix0,cellsHfixData) ;
std::vector<int> cellsHwel = indexC(cellsHwel0,cellsHwelData) ;
std::vector<int> cellsHdrn = indexC(cellsHdrn0,cellsHdrnData) ;
std::vector<int> cellsHghb = indexC(cellsHghb0,cellsHghbData) ;
std::vector<int> cellsHriv = indexC(cellsHriv0,cellsHrivData) ;
std::vector<int> cellsHrch = indexC(cellsHrch0,cellsHrchData) ;

//--------------------------- neural network ------------------------------------
scalar activateNNtransport(transportProperties.lookupOrDefault<scalar>("activateNNtransport",0));
scalar activateNNchemistry(transportProperties.lookupOrDefault<scalar>("activateNNchemistry",0));

//my_NN Cwgnn;
scalar nn_epoc(transportProperties.lookupOrDefault<scalar>("nn_epoc",6));
scalar nn_batch(transportProperties.lookupOrDefault<scalar>("nn_batch",16));
scalar nn_lr(transportProperties.lookupOrDefault<scalar>("nn_lr",1e-3));
scalar nn_minr(transportProperties.lookupOrDefault<scalar>("nn_minr",1e-3));
