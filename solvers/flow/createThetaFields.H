// base field

// get transport properties
Info << "Reading unsat properties\n" << endl;
//scalar n_vg(transportProperties.lookupOrDefault<scalar>("n_vg",1.5));
//dimensionedScalar alpha_vg(transportProperties.lookupOrDefault("alpha_vg",dimensionedScalar("alpha_vg",dimless/dimLength,10.0)));
volScalarField hp
(
    IOobject("hp",runTime.timeName(),mesh,IOobject::READ_IF_PRESENT,IOobject::AUTO_WRITE),mesh,dimensionedScalar("hp",dimLength,0)
);
volScalarField hp1=hp;
volScalarField hp_tmp=hp;
volScalarField sw_min
(
    IOobject("sw_min",runTime.constant(),mesh,IOobject::READ_IF_PRESENT,IOobject::AUTO_WRITE),mesh,dimensionedScalar("sw_min",dimless,1e-6)
);
volScalarField n_vg
(
    IOobject("n_vg",runTime.constant(),mesh,IOobject::READ_IF_PRESENT,IOobject::AUTO_WRITE),mesh,dimensionedScalar("n_vg",dimless,1.5)
);
volScalarField alpha_vg
(
    IOobject("alpha_vg",runTime.constant(),mesh,IOobject::READ_IF_PRESENT,IOobject::AUTO_WRITE),mesh,dimensionedScalar("alpha_vg",dimless/dimLength,1.)
);
volScalarField sw_max
(
    IOobject("sw_max",runTime.constant(),mesh,IOobject::READ_IF_PRESENT,IOobject::AUTO_WRITE),mesh,dimensionedScalar("sw_max",dimless,1-1e-4)
);
volScalarField m_vg = 1 - 1/n_vg;
volScalarField thetag ("thetag",unity);
volScalarField thetaw ("thetaw",unity);
volScalarField krw ("krw",unity);
volScalarField dkrwdS ("dkrwdS", unity);
volScalarField Se ("Se",unity);
volScalarField Se_tmp ("Se_tmp",unity);

surfaceScalarField krwf ("krwf",fvc::interpolate(krw,"krw"));
surfaceScalarField Lf = rhowf*Kf*krwf/muwf;
surfaceScalarField phiGr ("phiGr",(Lf * g) & mesh.Sf() );
//surfaceScalarField phiPcNull ("phiPc",0*phiGr);
label iterPicard=0;
scalar resPicard = 1000.;
Info << "end general fields theta, flowType" << flowType<<endl;

if (flowType == 3) // unsaturated
{
	//create head and calculate capillary pressure
	//hp.dimensions().reset(dimless);
	Se = pow(1.0+pow(-min(hp*alpha_vg,-1e-10),n_vg),-m_vg)-SMALL;
	//Se = neg(hp) *pow(1.0+pow(alpha_vg*mag(hp),n_vg),-m_vg) + pos(hp+VSMALL)*(1-SMALL);
	//hp.dimensions().reset(dimLength);
	sw = sw_min + (sw_max-sw_min)*Se;
	Info<<"n_vg "<<min(n_vg)<<" alpha "<<min(alpha_vg)<<" min Se "<<min(Se)<<" max Se "<<max(Se)<<" swmin 5 "<<sw_min[5]<<endl;
//Info << "end Se " << endl;

	//- Update effective saturation and capillary capacity
	// HERE Ch is called stor for compatibility
	stor.dimensions().reset(dimless/dimLength);
	//Ch == alpha*m*(Smax-Smin)/(1.0-m_vg)*pow(Se,1.0/m_vg)*pow(1.0-pow(Se,1.0/m_vg),m_vg) ;
	stor = eps*alpha_vg*m_vg*(sw_max-sw_min)/(1-m_vg)*pow(Se,1/m_vg)*pow(1-pow(Se,1/m_vg),m_vg);
	//stor.dimensions().reset(dimless/dimLength);
//Info << "end stor " << endl;

//- relative permeability (kr)
	krw = 1 * pow(Se,0.5) * pow(1-pow(1-pow(Se,1/m_vg),m_vg),2); //krbmax = 1
	krwf = fvc::interpolate(krw,"krw");
	//dkrwdS = 0.5 * (1-pow((1-pow(Se,1/m_vg)),m_vg));
	//dkrwdS *= ( 4 * pow(Se,1/m_vg-1/2) * pow( (1-pow(Se,1/m_vg)) , m_vg-1)) - ( pow((1-pow(Se,1/m_vg)),m_vg) -1) / pow(Se,0.5);
	//dkrwdS *= 1/(1 - sw_min);
	
	//- mobility and fractional flow 
	Mf = rhowf*mag(g)*Kf*krwf/muwf;// Mf alread exists
	Lf = rhowf*Kf*krwf/muwf;
	//- fluxes depending on saturation
	phiGr = (Lf * g) & mesh.Sf();

	Info << "end createThetaF \n" << endl;
}
