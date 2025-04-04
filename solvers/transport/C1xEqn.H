{
    // 
    sw = min(sw,0.999);
    sg = 1-sw;
	thetaw = eps*sw; 
	thetag = eps*sg;
if (diffusionEqn==0) {Deffg = vDg0 * tensor::I;}
if (diffusionEqn==1) {Deffg = (vDg0*pow(thetag,3.333)+Dw0*pow(thetaw,3.333))/pow(eps,2.) * tensor::I;}
if (diffusionEqn==2) {Deffg = (vDg0*pow(thetag,2.333)+Dw0*pow(thetaw,2.333))/pow(eps,2.) * tensor::I;}
if (diffusionEqn==3) {Deffg = vDg0 * thetag * tensor::I;}

// ---------------------GAS EQUATION ------------------------
// --------------------- internal conditions ------------------------
float time = mesh.time().value(); //const Time&
CSu=CSu*0; CSp=CSp*0; a1=1; double conc;//previous cells with values must be destroyed
size_t i;
if (cellsCgwel.size()>0)
{
	//get the concentration well data
	if (time>tCgwel) {
		output = getCbuffer("cgwel",itCgwel,cellsCgwel.size()); 
		tCgwel = output.t; cellsCgwelData=output.d; itCgwel +=1;
		tnext = min(tnext,tCgwel);
		}
	for (i=0; i<cellsCgwel.size(); i++) 
		{ 
		i1=cellsCgwel[i];v=mesh.V()[i1];
		if (qPwel[i]<0) {CSp[i1] = -qPwel[i]/v;} //in Sp term as discharge prop to local concentration, a positive Sp will make solution implicit
		else {CSu[i1] = -qPwel[i]/v*cellsCgwelData[i*4+2];} //the injected mass is the water flow x the conc
		}
}

//massTransfer_g = kmg * thetag * (Cg - Kdg*Cw);
massTransfer_g = kmg * (thetag*C1g+thetaw*C1w)/(1.0+thetag/thetaw*Kdg) - kmg * thetaw *C1w;
massTransfer_s = min(sw*eps,1.0-eps) * kms * (C1s - Kds*C1w);
thet1 = (1+thetag/thetaw*Kdg);
C1g.storePrevIter();
// ------------------- the equation ---------------------
fvScalarMatrix CgEqn
	(
		eps * sg * fvm::ddt(C1g) + eps * C1g * fvc::ddt(sg)
		+ fvm::div(phig, C1g, "div(phig,Cg)") // used only in multiphase
		- fvm::laplacian(Deffg, C1g, "laplacian(Deffg,Cg)")
		==
		- kmg * thetag/thet1*C1g - kmg*thetaw/thet1*C1w + kmg * thetaw *C1w
		- eps * sg * fvm::Sp(lbdag,C1g)
		- CSu - fvm::SuSp(CSp, C1g)
	);
//for cells of fixed conc (it is supposed that we don't put cfix at otuflow)
if (cellsCgfix.size()>0)
	{
	if (time>tCgfix) {
		output = getCbuffer("cgfix",itCgfix,cellsCgfix.size()); 
		tCgfix = output.t; cellsCgfixData=output.d; itCgfix +=1;
		tnext = min(tnext,tCgfix);
		}
	for (i=0; i<cellsCgfix.size(); i++) {CgEqn.setReference(cellsCgfix[i],cellsCgfixData[i*4+2]);}
	}

CgEqn.solve();

//Info<<"sw 20 "<<sw[20]<< " Deffg 10 "<<Deffg[10]<<" Deffg 20 "<<Deffg[20]<<" Cg 20 "<<Cg[20]<<endl;

//------------------- WATER EQUATION  -------------------------------
//Deffw = 1e-10;
//------------------ conditions ------------------
// for fixed head cells
	if (cellsHfix.size()>0)
	{
		if (time>tCwhfix) {
			output = getCbuffer("cwhfix",itCwhfix,cellsCwhfix.size()); 
			tCwhfix = output.t; cellsCwhfixData=output.d; itCwhfix +=1;
			tnext = min(tnext,tCwhfix);
			}
		for (size_t i=0; i<cellsHfix.size(); i++) 
			{ 
			i1=cellsHfix[i];v=mesh.V()[i1];//
			if (qHfix[i]<0) {CSp[i1] = -qHfix[i]/v;} // outflow discharge prop to local concentration
			else {CSu[i1] = -qHfix[i]/v*cellsCwhfixData[i*4+2];}
			}
	} Info<<"Chfix set"<<endl;

if (cellsCwwel.size()>0)
{
	//get the concentration well data
	if (time>tCwwel) {
		output = getCbuffer("cwwel",itCwwel,cellsCwwel.size()); 
		tCwwel = output.t; cellsCwwelData=output.d; itCwwel +=1;
		tnext = min(tnext,tCwwel);
		}
	for (i=0; i<cellsCwwel.size(); i++) 
		{ 
		i1=cellsCwwel[i];v=mesh.V()[i1];
		if (qHwel[i]<0) {CSp[i1] = -qHwel[i]/v;} //in Sp term as discharge prop to local concentration, a positive Sp will make solution implicit
		else {CSu[i1] = -qHwel[i]/v*cellsCwwelData[i*4+2];} //the injected mass is the water flow x the conc
		}
}

fvScalarMatrix CwEqn
	(
		eps * sw * fvm::ddt(C1w) + eps * C1w * fvc::ddt(sw)
		+ fvm::div(phiw, C1w, "div(phiw,Cw)")
		- fvm::laplacian(Deffw, C1w, "laplacian(Deffw,Cw)")
		==
		kmg * thetag/thet1*C1g.prevIter() + kmg*thetaw/thet1*C1w - kmg * thetaw *C1w 
		+ massTransfer_s
		- eps * sw * fvm::Sp(lbdaw,C1w)
		- CSu - fvm::SuSp(CSp, C1w)
	);
//for cells of fixed conc (it is supposed that we don't put cfix at otuflow)
if (cellsCwfix.size()>0)
	{
	if (time>tCwfix) {
		output = getCbuffer("cwfix",itCwfix,cellsCwfix.size()); 
		tCwfix = output.t; cellsCwfixData=output.d; itCwfix +=1;
		tnext = min(tnext,tCwfix);
		}
	for (i=0; i<cellsCwfix.size(); i++) {CwEqn.setReference(cellsCwfix[i],cellsCwfixData[i*4+2]);}
	}

CwEqn.solve();
//Cw = 1/Kds*Cg;

//----------------------------- solid  equation (no conditions) -------------------
fvScalarMatrix CsEqn
	(
		(1-eps) * fvm::ddt(C1s)
		 ==
		- massTransfer_s
		- (1-eps) * fvm::Sp(lbdas,C1s)
	);
CsEqn.solve();
Info<<"Cs solved "<<endl;
}
