Info<<"start 2phasefields"<<endl;
/*
autoPtr<incompressiblePhase> phaseg = incompressiblePhase::New(mesh,transportProperties,"g");
	volVectorField& Ug = phaseg->U();
	const dimensionedScalar& rhog = phaseg->rho();
	const dimensionedScalar& mug = phaseg->mu();
	phaseg->phi().writeOpt()=IOobject::NO_WRITE;
	Ug.readOpt()=IOobject::READ_IF_PRESENT;*/
	//Ug = Uw*0;//phig = phiw*0; //Ug is still needed to be read
dimensionedScalar rhog(transportProperties.lookupOrDefault("rhog",dimensionedScalar("rhog",dimensionSet(1,-3,0,0,0,0,0),1.2)));
dimensionedScalar mug(transportProperties.lookupOrDefault("mug",dimensionedScalar("mug",dimensionSet(1,-1,-1,0,0,0,0),1e-5)));
volVectorField Ug
(
    IOobject("Ug",runTime.timeName(),mesh,IOobject::READ_IF_PRESENT,IOobject::AUTO_WRITE),mesh,dimensionedVector("Ug",dimLength/dimTime,vector::zero)
);
surfaceScalarField phig ("phig",phiw*0);
//Info<<"sw "<<sw<<endl;
volScalarField p1=p;
volScalarField pw=p;
volScalarField p_tmp=p;
volScalarField pc0 = mag(g)*rhow/alpha_vg; //Info<<"pc0 "<<pc0<<endl;// alphavg 1/m -> Pa
Info<<"Sw min"<<sw_min<<" Sw max "<<sw_max<<endl;
// be careful, if g=0 pc0=0
Se = max(SMALL,min((sw-sw_min)/(sw_max-sw_min),1.-SMALL));
/*
if ((max(sw).value()-min(sw).value()>0.)||(pc0[0]==0.)) {Se = max(0.,min((sw-sw_min)/(sw_max-sw_min),1.));}
else { 
	Info<<"pc0 n m "<<pc0<<" "<<n_vg<<" "<<m_vg;
	Se =  pow(1.0+pow(max(-p/pc0,1.e-16),n_vg),-m_vg);Info<<"Se "<<Se[0]<<endl;
	//Se =  neg(hp) *pow(1.0+pow(alpha_vg*mag(hp),n_vg),-m_vg) + pos(hp+VSMALL)*(1-SMALL);
	sw = sw_min + (sw_max-sw_min)*Se;
	}
*/
volScalarField comp_g(IOobject("comp_g",runTime.constant(),mesh,IOobject::READ_IF_PRESENT,IOobject::AUTO_WRITE),mesh,dimensionedScalar("comp_g",dimensionSet(-1,1,2,0,0,0,0),1./atmPa));

volScalarField krg ("krg", pow(1-Se,0.5) * pow(1-pow(Se,1/m_vg),2*m_vg) );//Info<<"krg "<<krg<<endl;
volScalarField Seg = Se;
volScalarField sg (IOobject("sg",runTime.timeName(),mesh,IOobject::READ_IF_PRESENT,IOobject::AUTO_WRITE),mesh,dimensionedScalar("sg",dimless,0));
surfaceScalarField krgf ("krgf",fvc::interpolate(krg,"krg"));//Info<<"krgf "<<krgf<<endl;
volScalarField pc(IOobject("pc",runTime.constant(),mesh,IOobject::READ_IF_PRESENT,IOobject::AUTO_WRITE),mesh,dimensionedScalar("pc",dimensionSet(1,-1,-2,0,0,0,0),1.));
pc==pc0*pow(pow(Se+VSMALL,-1./m_vg)-1,1./n_vg);//Info<<"pc "<<pc<<endl;
volScalarField dpcdS ("dpcdS", pc0*0);//Info<<"dpcdS "<<dpcdS<<endl;
surfaceScalarField dpcdSf ("dpcdSf", fvc::interpolate(dpcdS) );
//volScalarField dkrgdS ("dkrgdS", unity);

//dpcdS.dimensions().reset(dimensionSet(1,-1,-2,0,0,0,0));

Info<<"end phase g"<<endl;
//- mobility computation 
surfaceScalarField Mwf ("Mwf",Kf*krwf/muwf);
surfaceScalarField Lwf ("Lwf",Mwf*rhowf);
surfaceScalarField Mgf ("Mgf",Kf*krgf/mug);
surfaceScalarField Lgf ("Lgf",Mgf*rhog);
surfaceScalarField Fwf ("Fwf",(krwf/muwf)/(krgf/mug + krwf/muwf));
volVectorField U = Ug + Uw;

volVectorField gradpc ("gradpc", fvc::grad(pc,"pc") );
surfaceVectorField gradpcf ("gradpcf", fvc::interpolate(gradpc,"gradpc") );
surfaceScalarField phi = phiw;
surfaceScalarField phiPc = ("phiPc",phiw*0);
surfaceScalarField phiP ("phiP", phiw);

if (flowType == 4) // 2 phase
	{
	// phase
	b.dimensions().reset(dimensionSet(-1,1,2,0,0,0,0));
	Mf.dimensions().reset(dimensionSet(-1,3,1,0,0,0,0)); // 
	Mf = Mwf + Kf*krgf/mug;//Info<<"Mf 0 1 "<<Mf[0]<<" "<<Mf[1]<<endl;
	Lf.dimensions().reset(dimensionSet(0,0,1,0,0,0,0));
	Lf = Lwf + Kf*rhog*krgf/mug;//Info<<"Lf "<<Lf[10]<<endl;
	surfaceScalarField phiPc ("phiPc", (Mwf * gradpcf) & mesh.Sf() );
	surfaceScalarField phiGr ("phiGr", (Lf * g) & mesh.Sf());//Info<<"phiGr 0 1 "<<phiGr[0]<<" "<<phiGr[1]<<endl;
	//recalculate fluxes (phiw, phig used for transport)
	//phi = phiGr - ((Mf * gradpc ) & mesh.Sf()) + phiPc*activateCapillarity;
    //phiw == Fwf*phi + (((Lwf/Lf)*g) & mesh.Sf()) + phiPc*activateCapillarity; // == to transfer also BC
	//phiw = phi*0;// (Lwf/Lf)*phiGr - ((Mwf * gradpc ) & mesh.Sf()) + phiPc*activateCapillarity;
Info<<"units phiGr"<<phiGr.dimensions()<<"units gradpc"<<gradpc.dimensions()<<"units phiPc"<<phiPc.dimensions()<<endl;
	Info<<" p ";for (i=0;i<8;i++) {Info<<p[i]<<" ";};Info<<"\n";
	Info<<" sw ";for (i=0;i<8;i++) {Info<<sw[i]<<" ";};Info<<"\n";
	Info<<" 1-Se ";for (i=0;i<8;i++) {Info<<1-Se[i]<<" ";};Info<<"\n";
	//Info<<" pc ";for (i=0;i<8;i++) {Info<<pc[i]<<" ";};Info<<"\n";
	Info<<" pc "<<pc;
	Info<<" phiGr ";for (i=0;i<8;i++) {Info<<phiGr[i]<<" ";};Info<<"\n";
	Info<<" Mwf ";for (i=0;i<8;i++) {Info<<Mwf[i]<<" ";};Info<<"\n";
	Info<<" grdpc ";for (i=0;i<8;i++) {Info<<gradpc[i]<<" ";};Info<<"\n";
	Info<<" phipc ";for (i=0;i<8;i++) {Info<<phiPc[i]<<" ";};Info<<"\n";
	}
Info<<"end Mf"<<endl;
//- compute fluxes

// Su
volScalarField::Internal pSu
(
   IOobject("pSu",runTime.timeName(),mesh,IOobject::NO_READ,IOobject::NO_WRITE),
		mesh,
		dimensionedScalar("zero",dimensionSet(0,0,-1,0,0,0,0),0.0) // the overall eq in in m3.s-1
);
volScalarField::Internal fSu
(
   IOobject("fSu",runTime.timeName(),mesh,IOobject::NO_READ,IOobject::NO_WRITE),
		mesh,
		dimensionedScalar("zero",dimensionSet(0,0,-1,0,0,0,0),0.0) // the overall eq in in m3.s-1
);
// -------------- conditions for pressure -------------
labelList cellsPfix0; std::vector<float> qPfix,cellsPfixData;
labelList cellsFfix0; std::vector<float> qFfix,cellsFfixData;
labelList cellsPwel0; std::vector<float> qPwel,cellsPwelData;
labelList cellsFwel0; std::vector<float> qFwel,cellsFwelData;
int itPfix,itPwel,itFfix,itFwel;float tPwel,tPfix,tFfix,tFwel;

std::vector<int>  cellsPfix,cellsPwel,cellsFfix,cellsFwel;
std::cout<<"start reader \n";
if (fexists(cur_dir+"/constant/options/pfix")) { 
std::cout<<"in pfix reader \n";
cellSet selectedPfix(mesh, "pfix");cellsPfix0 = selectedPfix.toc();qPfix.resize(cellsPfix0.size());
while (time>=tPfix) {output = getCbuffer("pfix",itPfix,cellsPfix0.size());itPfix +=1;tPfix = output.t;cellsPfixData=output.d;}
cellsPfix = indexC(cellsPfix0,cellsPfixData) ;
}
std::cout<<"end pfix reader \n";
if (fexists(cur_dir+"/constant/options/ffix")) { 
cellSet selectedFfix(mesh, "ffix");cellsFfix0 = selectedFfix.toc();qFfix.resize(cellsFfix0.size());
while (time>=tFfix) {output = getCbuffer("ffix",itFfix,cellsFfix0.size());itFfix +=1;tFfix = output.t;cellsFfixData=output.d;}
cellsFfix = indexC(cellsFfix0,cellsFfixData) ;
}
if (fexists(cur_dir+"/constant/options/pwel")) { 
cellSet selectedPwel(mesh, "pwel");cellsPwel0 = selectedPwel.toc();qPwel.resize(cellsPwel0.size());
while (time>=tPwel) {output = getCbuffer("pwel",itPwel,cellsPwel0.size());itPwel +=1;tPwel = output.t;cellsPwelData=output.d;}
cellsPwel = indexC(cellsPwel0,cellsPwelData) ;
}
if (fexists(cur_dir+"/constant/options/fwel")) { 
cellSet selectedFwel(mesh, "fwel");cellsFwel0 = selectedFwel.toc();qFwel.resize(cellsFwel0.size());
while (time>=tFwel) {output = getCbuffer("fwel",itFwel,cellsFwel0.size());itFwel +=1;tFwel = output.t;cellsFwelData=output.d;}
cellsFwel = indexC(cellsFwel0,cellsFwelData) ;
}
std::vector<int> bcPS(ncell,0);

