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
