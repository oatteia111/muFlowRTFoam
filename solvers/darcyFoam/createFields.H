
Info << nl << "Reading field p\n" << endl;
volScalarField h
(
    IOobject
    ("h",runTime.timeName(),mesh,IOobject::READ_IF_PRESENT,IOobject::AUTO_WRITE),
    mesh
);

volScalarField eps
(
    IOobject("eps",runTime.constant(),mesh,IOobject::READ_IF_PRESENT,IOobject::AUTO_WRITE),mesh,dimensionedScalar("eps",dimless,0.3)
);

Info << "Reading transportProperties\n" << endl;
IOdictionary transportProperties
(
    IOobject
    ("transportProperties",runTime.constant(), mesh,IOobject::MUST_READ_IF_MODIFIED,IOobject::NO_WRITE)
);

autoPtr<incompressiblePhase> phasew = incompressiblePhase::New(mesh,transportProperties,"w");
volVectorField& Uw = phasew->U();
const dimensionedScalar& rho = phasew->rho();
const dimensionedScalar& mu = phasew->mu();
phasew->phi().writeOpt()=IOobject::NO_WRITE;

Info<< "end read phase\n" << endl;

surfaceScalarField phiw
(
    IOobject
    ("phiw", runTime.timeName(), mesh,IOobject::READ_IF_PRESENT,IOobject::AUTO_WRITE),
	linearInterpolate(Uw) & mesh.Sf()
);

Info << "Reading field K\n" << endl;
volScalarField K
(
    IOobject
    ("K",runTime.constant(),mesh,IOobject::MUST_READ,IOobject::AUTO_WRITE),
    mesh
);

surfaceScalarField Kf(fvc::interpolate(K,"K"));
surfaceScalarField Mf ("Mf",rho*mag(g)*Kf/mu);
surfaceScalarField phiG("phiG",(rho * Kf/mu * g) & mesh.Sf());
