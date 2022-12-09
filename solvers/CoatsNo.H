// --- Coats (2003) IMPES STABILITY: CFL LIMIT
// - inertial part
{
    volScalarField dfw (
        "dfw",
        (dkrwdS*krg - dkrgdS*krw) /(Mmu*pow(krg,2)+2*krg*krw+1/Mmu*pow(krw,2))
    );
    dimensionedScalar smallRate("smallRate",dimVolume/dimTime, SMALL);
	//Info<<"sw[98] "<<sw[98]<<"sw[99] "<<sw[99]<<endl;
	//Info<<"dfw[98] "<<dfw[98]<<"dfw[99] "<<dfw[99]<<endl;
	//Info<<"krg[98] "<<krg[98]<<"krg[99] "<<krg[99]<<endl;
	//Info<<"krw[98] "<<krw[98]<<"krw[99] "<<krw[99]<<endl;
	//Info<<"dkrwdS[98] "<<dkrwdS[98]<<"dkrwdS[99] "<<dkrwdS[99]<<endl;
	//Info<<"dkrgdS[98] "<<dkrgdS[98]<<"dkrgdS[99] "<<dkrgdS[99]<<endl;
    // - gravity contribution (cahnged - to + from horgue because cfl became<0)
    dfw += (rhog-rhow)*fvc::surfaceSum(Kf*mag(mesh.Sf() & g))/fvc::surfaceSum(mag(phi)+smallRate)* (pow(krg,2)*dkrwdS/mug + pow(krw,2)*dkrgdS/muw)/(Mmu*pow(krg,2)+2*krg*krw+1/Mmu*pow(krw,2));

    volScalarField CFLCoats((runTime.deltaT()/eps)*dfw*fvc::surfaceSum(mag(phi)));
	//Info<<"dfw[98] "<<dfw[98]<<" dfw[99] "<<dfw[99]<<" CFLCoats[99] "<<CFLCoats[99]<<endl;

    // - capillarity part of CFL
    if(activateCapillarity)
    {
        CFLCoats += (runTime.deltaT()/eps)*2*mag(dpcdS)*(krg*krw/(muw*krg+mug*krw))*fvc::surfaceSum(Kf*mesh.magSf()/mag(mesh.delta()));
        // CFLCoats.ref() /= mesh.V();
    }

    CFLCoats.ref() /= mesh.V();
    CFLUse = gMax(cmptMax(CFLCoats)); //CFLCoats[0].component(0);
	Info<<"CFLCoats[6]"<<CFLCoats[6]<<" CFLuse "<<CFLUse<<" maxCo "<<maxCo<<endl;
    maxDeltaTFact = maxCo/(CFLUse + SMALL);
	Info<<"coats maxdTfact "<<maxDeltaTFact<<endl;

}
