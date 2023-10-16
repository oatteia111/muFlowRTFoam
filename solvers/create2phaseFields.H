Info<<"start 2phasefields"<<endl;
	autoPtr<incompressiblePhase> phaseg = incompressiblePhase::New(mesh,transportProperties,"g");
	volVectorField& Ug = phaseg->U();
	const dimensionedScalar& rhog = phaseg->rho();
	const dimensionedScalar& mug = phaseg->mu();
	phaseg->phi().writeOpt()=IOobject::NO_WRITE;
	Ug.readOpt()=IOobject::READ_IF_PRESENT;
	//Ug = Uw*0;//phig = phiw*0; //Ug is still needed to be read
	surfaceScalarField phig ("phig",phiw*0);
//Info<<"sw "<<sw<<endl;
Se = max(1e-6,min((sw-sw_min)/(sw_max-sw_min),1-1e-6));//Info<<"Se min"<<min(Se)<<" Se max "<<max(Se)<<endl;
volScalarField krg ("krg", pow(1-Se,0.5) * pow(1-pow(Se,1/m_vg),2*m_vg) );//Info<<"krg "<<krg<<endl;
surfaceScalarField krgf ("krgf",fvc::interpolate(krg,"krg"));//Info<<"krgf "<<krgf<<endl;
//dimensionedScalar pc0 ("pc0",dimensionSet(1,-1,-2,0,0,0,0));
volScalarField pc0 = mag(g)*rhow/alpha_vg; //Info<<"pc0 "<<pc0<<endl;// alphavg 1/m -> Pa
volScalarField dpcdS ("dpcdS", - 1/(n_vg*m_vg) * pc0/(sw_max - sw_min) * pow(pow(Se,-1/m_vg)-1,(1/n_vg)-1) * pow(Se,-(1+m_vg)/m_vg) );//Info<<"dpcdS "<<dpcdS<<endl;
volScalarField dkrgdS ("dkrgdS", unity);

	Info << "pc0 dim "<<pc0.dimensions()<<" dpcdS dim" << dpcdS.dimensions() << endl;
	dpcdS.dimensions().reset(dimensionSet(1,-1,-2,0,0,0,0));

	dimensionedScalar Mmu (muw/mug);//viscosity ratio
	surfaceScalarField Fwf ("Fwf",(krwf/muw) / ( (krgf/mug) + (krwf/muw) ) );

Info<<"end phase g"<<endl;
//- mobility computation 
surfaceScalarField Mwf ("Mwf",Kf*krwf/muw);
surfaceScalarField Lwf ("Lwf",Mwf*rhow);
if (flowType == 4) // 2 phase
	{
	// phase
	U = Ug + Uw;

	b.dimensions().reset(dimensionSet(-1,1,2,0,0,0,0));
	Mf.dimensions().reset(dimensionSet(-1,3,1,0,0,0,0)); // 
	Mf = Mwf + Kf*krgf/mug;Info<<"Mf 0 1 "<<Mf[0]<<" "<<Mf[1]<<endl;
	Lf.dimensions().reset(dimensionSet(0,0,1,0,0,0,0));
	Lf = Lwf + Kf*rhog*krgf/mug;//Info<<"Lf "<<Lf[10]<<endl;
	}
Info<<"end Mf"<<endl;
//- compute fluxes
surfaceVectorField gradpc ("gradpc", fvc::interpolate(dpcdS*fvc::grad(sw),"pc") );
surfaceScalarField phiPc ("phiPc", (Mwf * gradpc) & mesh.Sf() );
phiGr = (Lf * g) & mesh.Sf();Info<<"phiGr 0 1 "<<phiGr[0]<<" "<<phiGr[1]<<endl;
surfaceScalarField phi ("phi", phiw);
Info<<"units phiGr"<<phiGr.dimensions()<<"units gradpc"<<gradpc.dimensions()<<"units phiPc"<<phiPc.dimensions()<<endl;

// Su
volScalarField::Internal pSu
(
   IOobject("pSu",runTime.timeName(),mesh,IOobject::NO_READ,IOobject::NO_WRITE),
		mesh,
		dimensionedScalar("zero",dimensionSet(0,0,-1,0,0,0,0),0.0) // the overall eq in in m3.s-1
);
// -------------- conditions for pressure -------------
labelList cellsPfix0; std::vector<float> qPfix;std::vector<float> condPfix;std::vector<float> cellsPfixData;
labelList cellsPwel0; std::vector<float> qPwel;std::vector<float> condPwel;std::vector<float> cellsPwelData;
int itPfix=0;int itPwel=0;float tPwel=0;float tPfix=0;

if (fexists(cur_dir+"/constant/options/pfix")) {cellSet selectedPfix(mesh, "pfix");cellsPfix0 = selectedPfix.toc();}
if (fexists(cur_dir+"/constant/options/pwel")) {cellSet selectedPwel(mesh, "pwel");cellsPwel0 = selectedPwel.toc();}

while (time>=tPfix) {output = getCbuffer("pfix",itPfix,cellsPfix0.size()); itPfix+=1;tPfix = output.t;cellsPfixData=output.d;}
while (time>=tPwel) {output = getCbuffer("pwel",itPwel,cellsPwel0.size()); itPwel+=1;tPwel = output.t;cellsPwelData=output.d;}

qPfix.resize(cellsPfix0.size());
qPwel.resize(cellsPwel0.size());
condPfix.resize(cellsPfix0.size());
condPwel.resize(cellsPwel0.size());

std::vector<int> cellsPfix = indexC(cellsPfix0,cellsPfixData) ;
std::vector<int> cellsPwel = indexC(cellsPwel0,cellsPwelData) ;

