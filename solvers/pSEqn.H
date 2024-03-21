// combine p and S equations
{

// -################# - conditions  ---####################---
	pSu=pSu*0;size_t i; //previous cells with values must be destroyed
	float time = mesh.time().value();float f; 
	icount = 0;	//get the ghb conditions
	if (cellsPwel.size()>0)
		{
		if (time>=tPwel) {
			output = getCbuffer("pwel",itPwel,cellsPwel.size()); 
			tPwel = output.t; cellsPwelData=output.d; itPwel +=1;
			tnext = tPwel;
			}
		for (i=0; i<cellsPwel.size(); i++) {condPwel[i] = static_cast<int>(cellsPwelData[i*4+2]);}
		}
	Info<<"Pfwel read ";

	fvScalarMatrix pEqn
        (
            fvm::laplacian(-Mf, p) + fvc::div(phiGr)
            + fvc::div(phiPc)*activateCapillarity
			- pSu
        );
		
if (cellsPfix.size()>0)
	{
	if (time>=tPfix) {
		output = getCbuffer("pfix",itPfix,cellsPfix.size()); 
		tPfix = output.t; cellsPfixData=output.d; itPfix +=1;
		tnext= tPfix;
		}
	for (i=0; i<cellsPfix.size(); i++) {pEqn.setReference(cellsPfix[i],static_cast<int>(cellsPfixData[i*4+2]));}
	}

	fvScalarMatrix SEqn
        (
            eps*fvm::ddt(sw) + fvc::div(phiw) 
            //==
			//b*fvOptions(sw,sw)
        );

	pEqn.solve();
	//for (j=0; j<3;j++) {Info<<"p after peq "<<p[j]<<" sw "<<sw[j]<<endl;}

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