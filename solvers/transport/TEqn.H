
// diffusion+dispersion
dimensionedScalar USMALL("",dimLength/dimTime,1e-36);
//Deffw = (alphaT * mag(vp))*tensor::I + (alphaL-alphaT) * (vp*vp)/(mag(vp)+USMALL); // Ut rpelaced by vp
thetaw = eps*sw;
Deffw = (alphaT * mag(Uw) + Dw0* pow(thetaw,3.333) /pow(eps,2.))*tensor::I + (alphaL-alphaT) * (Uw*Uw)/(mag(Uw)+USMALL);//Info<<"Deffw "<<Deffw<<endl;

forAll(Deffw.boundaryField(),patchi)
	{ //- Dispersion equal to zero on patch where flux is fixed
	if (C.boundaryField().types()[patchi] == "fixedFlux") Deffw.boundaryFieldRef()[patchi] = tensor(zero());
	}
lbdaT = thetaw*lbdaTw + (1-eps)*lbdaTs;

float time = mesh.time().value(); //const Time&
// --------------------- internal conditions ------------------------
if (cellsHwel.size()>0)
	{
		//get the temp well data
		if (time>tTwel) {
			output = getCbuffer("twel",itTwel,cellsTwel.size()); 
			tTwel = output.t; cellsTwelData=output.d; itTwel +=1;
			tnext = min(tnext,tTwel);
			}
		for (size_t i=0; i<cellsTwel.size(); i++) 
			{ 
			i1=cellsTwel[i];v=mesh.V()[i1];
			if (qHwel[i]<0) {TSp[i1] = -rhow[i1]*qHwel[i]/v;} //in Sp term as discharge prop to local concentration, a positive Sp will make solution implicit
			else {TSu[i1] = -rhow[i1]* qHwel[i]/v*cellsTwelData[i*4+2];} //the injected mass is the water flow x the conc
			//Info<<"cellw "<<cellsTwel[i]<<" q "<<qHwel[i]<<" tloc "<<T[i1]<<" Su "<<TSu[i1]<<" Sp "<<TSp[i1]<<endl;
			}
	}

	// for fixed head cells
	if (cellsHfix.size()>0)
	{
		if (time>tThfix) {
			output = getCbuffer("thfix",itThfix,cellsThfix.size()); 
			tThfix = output.t; cellsThfixData=output.d; itThfix +=1;
			tnext = min(tnext,tThfix);
			}
		for (size_t i=0; i<cellsHfix.size(); i++) 
			{ 
			i1=cellsHfix[i];v=mesh.V()[i1];//
			if (qHfix[i]<0) {TSp[i1] = -rhow[i1]*qHfix[i]/v;} // outflow discharge prop to local concentration
			else {TSu[i1] = -rhow[i1]*qHfix[i]/v*cellsThfixData[i*4+2];}
			//Info<<"cell thfix "<<cellsThfix[i]<<" q "<<qHfix[i]<<" tloc "<<T[i1]<<" Su "<<TSu[i1]<<" Sp "<<TSp[i1]<<endl;
			}
	}

	// fixed cells done with ~ghb with very high conductance
	if (cellsTfix.size()>0)
	{
		if (time>tTfix) {
			output = getCbuffer("tfix",itTfix,cellsTfix.size()); 
			tTfix = output.t; cellsTfixData=output.d; itTfix +=1;
			tnext = min(tnext,tTfix);
			}
		for (size_t i=0; i<cellsTfix.size(); i++) 
			{ 
			i1=cellsTfix[i];v=mesh.V()[i1];TSu[i1] = -cellsTfixData[i*4+2]*1e5; TSp[i1] = 1e5;// high C conductance
			}
	}

T.storePrevIter();
surfaceScalarField rhophi = rhowf* phiw;

Info << "before Teq "<<endl;

fvScalarMatrix TEqn
	(
		rhow*cpw * thetaw * fvm::ddt(T) + rhow*cpw * T * fvc::ddt(thetaw) + rhos*cps * (1-eps) * fvm::ddt(T) 
		+ fvc::div(cpw * rhophi, T, "div(phiw,T)") 
		- fvm::laplacian(eps*rhow*cpw * Deffw, T, "laplacian(Deffw,T)")
		- fvm::laplacian(lbdaT, T, "laplacian(lbdaT,T)")
		+ cpw*fvm::SuSp(TSp, T)
		== - cpw*TSu 
	);
	
if (inactiveCells.size()+dryCells.size()>0) 
{
	Info<<"in inact/dry Cells, nb inact "<<inactiveCells.size()<<" dry size "<<dryCells.size()<<endl;
	labelHashSet inactiveCellsSet(inactiveCells);	
	forAll(mesh.owner(),iface)
	{
		if ((inactiveCellsSet.found(TEqn.lduAddr().lowerAddr()[iface]))||(inactiveCellsSet.found(TEqn.lduAddr().upperAddr()[iface])))
		{
			TEqn.upper()[iface] = 0.0;TEqn.lower()[iface] = 0.0;
		}
	}
	for (label cellI : inactiveCells) {TEqn.diag()[cellI] = 1.0;TEqn.source()[cellI] = 0;}
	
	labelList drfList = dryFaces.shrink(); // neighb is low
	for (label faceI : drfList) {TEqn.lower()[faceI] = 0.0;TEqn.upper()[faceI] = 0.0;}
	labelList drcList = dryCells.shrink();
	for (label cellI : drcList) {TEqn.diag()[cellI] = 1.0;TEqn.source()[cellI] = 0;}
}

/*
//get the fixed temperature conditions
if (cellsTfix.size()>0)
{
	if (time>tTfix) {
		output = getCbuffer("tfix",itTfix,cellsTfix.size()); 
		tTfix = output.t; cellsTfixData=output.d; itTfix +=1;
		tnext = min(tnext,tTfix);
		}
	for (size_t i=0; i<cellsTfix.size(); i++) 
		{TEqn.setReference(cellsTfix[i],cellsTfixData[i*4+2]);} 
	Info<<"tfix in Eqn "<<time<<" val "<<cellsTfixData[2]<<endl;
}
*/

TEqn.solve();

//calculating dT
for (i=0;i<ncell;i++) {dT= max(dT,T()[i]-T.prevIter()[i]);}



