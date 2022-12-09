
// diffusion+dispersion
dimensionedScalar USMALL("",dimLength/dimTime,1e-36);
//Deffw = (alphaT * mag(vp))*tensor::I + (alphaL-alphaT) * (vp*vp)/(mag(vp)+USMALL); // Ut rpelaced by vp
thetaw = eps*sw;
Deffw = (alphaT * mag(Uw) + Dw0* pow(thetaw,3.333) /pow(eps,2.))*tensor::I + (alphaL-alphaT) * (Uw*Uw)/(mag(Uw)+USMALL);//Info<<"Deffw "<<Deffw<<endl;

forAll(Deffw.boundaryField(),patchi)
	{ //- Dispersion equal to zero on patch where flux is fixed
	if (Cw.boundaryField().types()[patchi] == "fixedFlux") Deffw.boundaryFieldRef()[patchi] = tensor(zero());
	}
Cw.storePrevIter();
//double sC=0;//
//std::vector<double> sC(ncell_lay,0);
//for (i=0;i<ncell_lay;i++) {sC+=abs(Cw()[i]);mxC=max(mxC,Cw()[i]);}
fvScalarMatrix CwEqn
	(
		R * thetaw * fvm::ddt(Cw) + R * Cw * fvc::ddt(thetaw) 
		+ fvm::div(phiw, Cw, "div(phiw,Cw)") 
		- fvm::laplacian(Deffw, Cw, "laplacian(Deffw,Cw)")
		+ thetaw * fvm::Sp(lbdaw,Cw)
		==
		fvOptions(sw,Cw) 
	);
fvOptions.constrain(CwEqn);
CwEqn.solve();
//Info <<"resid "<<CwEqn.residual()()[67171]<<endl;//.value()<<endl;
//Info << "Ceqn solve, alphaL "<<alphaL<<" theta 5 "<<thetaw[5]<<" Deffw 5"<<Deffw[5]<<endl;
//scalar dCdTmax = max(mag(fvc::ddt(Cw))).value();Info<<"dCdTmax "<<dCdTmax<<endl;
//scalar dC = dCdTmax*runTime.deltaTValue()*1;
//Info<<"max dCw " <<max(mag(Cw-Cw.prevIter()))<<endl;
double dC=0;double mxC=0;double sC=0;double sC1=0;
a *= 0;
for (i=0;i<ncell_lay;i++) {mxC=max(mxC,Cw()[i]);sC+=Cw()[i];}
Cw=max(Cw,a);
for (i=0;i<ncell_lay;i++) {dC= max(dC,Cw()[i]-Cw.prevIter()[i]);sC1+=Cw()[i];}
Cw *= sC/(sC1+SMALL); // removing negative value increase artifially the total
dC = dC/(mxC+SMALL);
Info<<"ncell "<<ncell_lay<<"dC "<<dC<<endl;
//scalar dC = sC1-sC;
//Info<<"max dCw "<<Cw[5]<<" "<<Cw.prevIter()[5]<<" "<<Cw[5]-Cw.prevIter()[5]<<endl;;//.internalField()


