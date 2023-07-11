// combine p and S equations
{

// -################# - conditions  ---####################---
	pSu*=0; //previous cells with values must be destroyed
	float time = mesh.time().value();float f; 
	icount = 0;	//get the ghb conditions
	if (cellsPw.size()>0)
	{
		float tnow = cellsPwData[2]*day;//time in files is in days
		while (time>=tnow) {icount += cellsPw.size(); tnow = cellsPwData[2+icount*4]*day; }
		if (icount>0) {icount -= cellsPw.size();};
		for (int i=0; i<cellsPw.size(); i++)
			{v=mesh.V()[cellsPw[i]];qPw[i] = cellsPwData[2+(icount+i)*4+2]/day; pSu[cellsPw[i]] = qPw[i]/v;}
	} Info<<"Pfw read ";

	fvScalarMatrix pEqn
        (
            fvm::laplacian(-Mf, p) + fvc::div(phiGr)
            + fvc::div(phiPc)*activateCapillarity
			- pSu
        );
		
	if (cellsPfx.size()>0)
	{
		float tnow = cellsPfxData[2]*day; 
		while (time>=tnow) {icount += cellsPfx.size(); tnow = cellsPfxData[2+icount*4]*day; }
		if (icount>0) {icount -= cellsPfx.size();} 
		for (int i=0; i<cellsPfx.size(); i++) {pEqn.setReference(cellsPfx[i],cellsPfxData[2+(icount+i)*4+2]);}
	}

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