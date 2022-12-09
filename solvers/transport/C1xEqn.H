{
    // 
    sw = min(sw,0.999);
    sg = 1-sw;
	thetaw = eps*sw; 
	thetag = eps*sg;
	Deffg = Dg0 * pow(thetag,3.33)/pow(eps,2.)*tensor::I;

	//massTransfer_g = kmg * thetag * (Cg - Kdg*Cw);
	massTransfer_g = kmg * (thetag*Cg+thetaw*Cw)/(1+thetag/thetaw*Kdg) - kmg * thetaw *Cw;
	massTransfer_s = min(sw*eps,1-eps) * kms * (Cs - Kds*Cw);
	thet1 = (1+thetag/thetaw*Kdg);
	Cg.storePrevIter();
	// gaz		
	fvScalarMatrix CgEqn
		(
			eps * sg * fvm::ddt(Cg) + eps * Cg * fvc::ddt(sg)
			+ fvm::div(phig, Cg, "div(phig,Cg)") // used only in multiphase
			- fvm::laplacian(Deffg, Cg, "laplacian(Deffg,Cg)")
			==
			- kmg * thetag/thet1*Cg - kmg*thetaw/thet1*Cw + kmg * thetaw *Cw
			- eps * sg * fvm::Sp(lbdag,Cg)
			+ fvOptions(sw,Cg)
		);
	CgEqn.solve();
	Info<<"sw 20 "<<sw[20]<< " Deffg 10 "<<Deffg[10]<<" Deffg 20 "<<Deffg[20]<<" Cg 20 "<<Cg[20]<<endl;

	//water
	//Deffw = 1e-10;
	fvScalarMatrix CwEqn
		(
			eps * sw * fvm::ddt(Cw) + eps * Cw * fvc::ddt(sw)
			+ fvm::div(phiw, Cw, "div(phiw,Cw)")
			- fvm::laplacian(Deffw, Cw, "laplacian(Deffw,Cw)")
			==
			kmg * thetag/thet1*Cg.prevIter() + kmg*thetaw/thet1*Cw - kmg * thetaw *Cw 
			+ massTransfer_s
			- eps * sw * fvm::Sp(lbdaw,Cw)
			+ fvOptions(sw,Cw)
		);
	CwEqn.solve();
	//Cw = 1/Kds*Cg;
	
   //solid
	fvScalarMatrix CsEqn
		(
			(1-eps) * fvm::ddt(Cs)
			 ==
			- massTransfer_s
			- (1-eps) * fvm::Sp(lbdas,Cs)
		);
	CsEqn.solve();
	Info<<"Cs solved "<<endl;
}
