
// diffusion+dispersion
dimensionedScalar USMALL("",dimLength/dimTime,1e-36);
//Deffw = (alphaT * mag(vp))*tensor::I + (alphaL-alphaT) * (vp*vp)/(mag(vp)+USMALL); // Ut rpelaced by vp
thetaw = eps*sw;
Deffw = (alphaT * mag(Uw) + Dw0* pow(thetaw,3.333) /pow(eps,2.))*tensor::I + (alphaL-alphaT) * (Uw*Uw)/(mag(Uw)+USMALL);//Info<<"Deffw "<<Deffw<<endl;

forAll(Deffw.boundaryField(),patchi)
	{ //- Dispersion equal to zero on patch where flux is fixed
	if (C.boundaryField().types()[patchi] == "fixedFlux") Deffw.boundaryFieldRef()[patchi] = tensor(zero());
	}
C.storePrevIter();
//double sC=0;//
//std::vector<double> sC(ncell_lay,0);
//for (i=0;i<ncell_lay;i++) {sC+=abs(Cw()[i]);mxC=max(mxC,Cw()[i]);}
fvScalarMatrix CEqn
	(
		R * thetaw * fvm::ddt(C) + R * C * fvc::ddt(thetaw) 
		+ fvm::div(phiw, C, "div(phiw,C)") 
		- fvm::laplacian(Deffw, C, "laplacian(Deffw,C)")
		+ thetaw * fvm::Sp(lbdaw,C)
		==
		fvOptions(sw,C) 
	);
fvOptions.constrain(CEqn);
CEqn.solve();
//Info <<"resid "<<CwEqn.residual()()[67171]<<endl;//.value()<<endl;
//Info << "Ceqn solve, alphaL "<<alphaL<<" theta 5 "<<thetaw[5]<<" Deffw 5"<<Deffw[5]<<endl;
//scalar dCdTmax = max(mag(fvc::ddt(Cw))).value();Info<<"dCdTmax "<<dCdTmax<<endl;
//scalar dC = dCdTmax*runTime.deltaTValue()*1;
//Info<<"max dCw " <<max(mag(Cw-Cw.prevIter()))<<endl;

//calculating dC, removing negative C and compensating mass
double dC=0;double mxC=0;double sC=0;double sC1=0;
a *= 0;
for (i=0;i<ncell;i++) {mxC=max(mxC,C()[i]);sC+=C()[i];}
C=max(C,a); // removing negative
for (i=0;i<ncell;i++) {dC= max(dC,C()[i]-C.prevIter()[i]);sC1+=C()[i];}
C *= sC/(sC1+SMALL); // removing negative value increase artifially the total
dC = dC/(mxC+SMALL);


