
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
	
/*forAll(mesh.cells(), celli)
    {
		const labelList& cFaces = mesh.cells()[celli];
        forAll(cFaces, cFacei){
			Info<<"cell "<<celli<<" faces "<<cFaces[cFacei]<<" phiw "<<phiw[cFaces[cFacei]]<<" "<<mesh.neighbour()[cFaces[cFacei]]<<endl;
			} 
	}*/
//this is the flux divergence
//forAll(dvPhi,i) {if (dvPhi()[i]>1e-3) {TSp[i] = -dvPhi()[i];Info<<"div pos at "<<i<<" dv "<<dvPhi()[i]<<endl;} } //if the div is too positive there is imposed discharge here

T.storePrevIter();
surfaceScalarField rhophi = rhowf* phiw;

Info << "before Teq "<<endl;

fvScalarMatrix TEqn
	(
		rhow*cpw * thetaw * fvm::ddt(T) + rhow*cpw * T * fvc::ddt(thetaw) + rhos*cps * (1-eps) * fvm::ddt(T) 
		+ fvm::div(cpw * rhophi, T, "div(phiw,T)") 
		- fvm::laplacian(eps*rhow*cpw * Deffw, T, "laplacian(Deffw,T)")
		- fvm::laplacian(lbdaT, T, "laplacian(lbdaT,T)")
		== - cpw*TSu - cpw*fvm::SuSp(TSp, T)
	);
	
//get the fixed temperature conditions
if (cellsTfix.size()>0)
{
	if (time>=tTfix) {
		output = getCbuffer("tfix",itTfix,cellsTfix.size()); 
		tTfix = output.t; cellsTfixData=output.d; itTfix +=1;
		tnext = min(tnext,tTfix);
		}
	for (size_t i=0; i<cellsTfix.size(); i++) 
		{TEqn.setReference(cellsTfix[i],cellsTfixData[i*4+2]);} 
	Info<<"tfix in Eqn "<<time<<" val "<<cellsTfixData[2]<<endl;
}

TEqn.solve();

//calculating dT
for (i=0;i<ncell;i++) {dT= max(dT,T()[i]-T.prevIter()[i]);}



