// combine p and S equations
{

	//for (j=0; j<3;j++) {Info<<"p before peq "<<p[j]<<" sw "<<sw[j]<<endl;}
	fvScalarMatrix pEqn
        (
            fvm::laplacian(-Mf, p) + fvc::div(phiGr)
            + fvc::div(phiPc)*activateCapillarity
            //==
			//b*fvOptions(sw,p)
        );
    //fvOptions.constrain(pEqn);
	pEqn.solve();
	//for (j=0; j<3;j++) {Info<<"p after peq "<<p[j]<<" sw "<<sw[j]<<endl;}

	fvScalarMatrix SEqn
        (
            eps*fvm::ddt(sw) + fvc::div(phiw) 
            //==
			//b*fvOptions(sw,sw)
        );
    SEqn.solve();
	//for (j=0; j<3;j++) {Info<<"p after seq "<<p[j]<<" sw "<<sw[j]<<endl;}
	
	#include "update2phaseFields.H"		
	//for (j=0; j<3;j++) {Info<<"p after upd "<<p[j]<<" sw "<<sw[j]<<endl;}

	//Info<<"p solved"<<endl;
    phi = pEqn.flux() + phiGr + phiPc*activateCapillarity;

    U = fvc::reconstruct(phi);
    U.correctBoundaryConditions();

    phiw == Fwf*phi + (((Lwf-Fwf*Lf)*g) & mesh.Sf()) + (1-Fwf)*phiPc*activateCapillarity;
    phig == phi - phiw;
	for (j=0; j<3;j++) {Info<<"phiw "<<phiw[j]<<endl;}

    Uw = fvc::reconstruct(phiw);
    Ug = U-Uw;

    Ug.correctBoundaryConditions();
    Uw.correctBoundaryConditions();  
}