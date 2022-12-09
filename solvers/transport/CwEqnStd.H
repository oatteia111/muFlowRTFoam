//*****************  steady flow *****************
simpleControl simple(mesh);  // new 11/6/21

// diffusion+dispersion
dimensionedScalar USMALL("",dimLength/dimTime,1e-36);
//Deffw = (alphaT * mag(vp))*tensor::I + (alphaL-alphaT) * (vp*vp)/(mag(vp)+USMALL); // Ut rpelaced by vp
Deffw = (alphaT * mag(Uw))*tensor::I + (alphaL-alphaT) * (Uw*Uw)/(mag(Uw)+USMALL);
forAll(Deffw.boundaryField(),patchi)
	{ //- Dispersion equal to zero on patch where flux is fixed
	if (Cw.boundaryField().types()[patchi] == "fixedFlux") Deffw.boundaryFieldRef()[patchi] = tensor(zero());
	}
scalarField resid = Cw*0;
//resid.dimensions().reset(dimensionSet(1,-3,-1,0,0,0,0));
int cnt=0;
//volScalarField qopt=sw*0;
//qopt.dimensions().reset(dimensionSet(0,0,-1,0,0,0,0));
while (simple.loop(runTime))
{
	while (simple.correctNonOrthogonal())
    {
	fvScalarMatrix CwEqnStd
		(
			fvm::div(phiw, Cw, "div(phiw,Cw)") 
			- fvm::laplacian(Deffw, Cw, "laplacian(Deffw,Cw)")
			+ eps * sw * fvm::Sp(lbdaw,Cw)
			//+ fvm::Sp(qopt,Cw)
			==
			fvOptions(sw,Cw) 
		);
	//Cw.storePrevIter();
	fvOptions.constrain(CwEqnStd);
	//Info<<"Cw time "<<runTime.timeName()<<" sum resid "<<sum(CwEqnStd.residual())<<endl;
	CwEqnStd.solve();
	//Info <<"resid "<<CwEqnStd.residual()()[67171]<<endl;//.value()<<endl;
	//fvOptions.correct(Cw);
	//Cw = Cw*0.2+Cw.prevIter()*0.8;
	Cw.relax();
	//Info << "Ceqn solve "<<endl;
	resid = CwEqnStd.residual();
	cnt = 1;
	}
}
//if (runTime.time().value()-cnt>10000) {
Cw.write();h.write();phiw.write();
std::ofstream cresFile(cur_dir/"Cwresiduals");
//for (i==0;i<cactive.size();i++) {outFile << CwEqnStd.residual()()[i] << "\n";}
for (const auto &x : resid) {cresFile << x << "\n";}
cresFile.close();
		//Info <<"resid "<<CwEqnStd.residual()().write()<<endl;

