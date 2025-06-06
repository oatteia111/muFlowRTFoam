
// diffusion+dispersion
dimensionedScalar USMALL("",dimLength/dimTime,1e-36);
//Deffw = (alphaT * mag(vp))*tensor::I + (alphaL-alphaT) * (vp*vp)/(mag(vp)+USMALL); // Ut rpelaced by vp
//if (flowType==2) {thetaw = sw;}
//else {thetaw = eps*sw;}
thetaw=eps*sw;
Deffw = (alphaT * mag(Uw) + Dw0* pow(thetaw,3.333) /pow(eps,2.))*tensor::I + (alphaL-alphaT) * (Uw*Uw)/(mag(Uw)+USMALL);//Info<<"Deffw "<<Deffw<<endl;
//forAll(Deffw(), i)  { if (sw[i]<sw_min[i]*1.5) {Deffw[i] *=1e2;} }; 

/*forAll(Deffw.boundaryField(),patchi)
	{ //- Dispersion equal to zero on patch where flux is fixed
	if (C.boundaryField().types()[patchi] == "fixedFlux") Deffw.boundaryFieldRef()[patchi] = tensor(zero());
	}*/
//double sC=0;//
//std::vector<double> sC(ncell_lay,0);
//for (i=0;i<ncell_lay;i++) {sC+=abs(Cw()[i]);mxC=max(mxC,Cw()[i]);}

// ---------------  global view including NN -------------
// Normally, if the conditions don't change and the model is trained we use NN, if not we use the eqs

// --------------------- internal conditions ------------------------
	float time = mesh.time().value(); //const Time&
	CSu=CSu*0; CSp=CSp*0; a1=1; double conc;//previous cells with values must be destroyed
	size_t i;
	if (cellsHwel.size()>0)
	{
		//get the concentration well data
		if (time>tCwel) {
			output = getCbuffer("cwel",itCwel,cellsCwel.size()); 
			tCwel = output.t; cellsCwelData=output.d; itCwel +=1;
			tnext = min(tnext,tCwel);
			}
		for (i=0; i<cellsCwel.size(); i++) 
			{ 
			i1=cellsCwel[i];v=mesh.V()[i1];
			if (qHwel[i]<0) {CSp[i1] = -qHwel[i]/v;} //in Sp term as discharge prop to local concentration, a positive Sp will make solution implicit
			else {CSu[i1] = -qHwel[i]/v*cellsCwelData[i*4+2];} //the injected mass is the water flow x the conc
			}
	}
	// Hghb because Cghb may not be set here
	if (cellsHghb.size()>0)
	{
		if (time>tCghb) {
			output = getCbuffer("cghb",itCghb,cellsCghb.size()); 
			tCghb = output.t; cellsCghbData=output.d; itCghb +=1;
			tnext = min(tnext,tCghb);
			}
		for (size_t i=0; i<cellsHghb.size(); i++) 
			{ 
			i1=cellsHghb[i];v=mesh.V()[i1];//
			if (qHghb[i]<0) {CSp[i1] = -qHghb[i]/v;} // outflow discharge prop to local concentration
			else {CSu[i1] = -qHghb[i]/v*cellsCghbData[i*4+2];}
			//Info<<"ghb "<< i<<" i1 " <<i1<<" q "<<qHghb[i]<<" conc "<<cellsCghbData[i*4+2]<<endl;
			}
	}
	//store the chfix conditions data
	// for fixed head cells
	if (cellsHfix.size()>0)
	{
		if (time>tChfix) {
			output = getCbuffer("chfix",itChfix,cellsChfix.size()); 
			tChfix = output.t; cellsChfixData=output.d; itChfix +=1;
			tnext = min(tnext,tChfix);
			}
		for (size_t i=0; i<cellsHfix.size(); i++) 
			{ 
			i1=cellsHfix[i];v=mesh.V()[i1];//
			if (qHfix[i]<0) {CSp[i1] = -qHfix[i]/v;} // outflow discharge prop to local concentration
			else {CSu[i1] = -qHfix[i]/v*cellsChfixData[i*4+2];}
			}
	} Info<<"Chfix set"<<endl;
	// for drain cells (only out), same as above, except that drains are constant in time
	if (cellsHdrn.size()>0)
	{
		for (size_t i=0; i<cellsHdrn.size(); i++) 
			{ 
			i1=cellsHdrn[i];v=mesh.V()[i1];//
			if (qHdrn[i]<0) {CSp[i1] = -qHdrn[i]/v;} // outflow qHdrn<0 discharge prop to local concentration
			}
	}

	//recharge
	if (cellsCrch.size()>0) // recharge is different, required only if C>0
	{
		if (time>tCrch) {
			output = getCbuffer("crch",itCrch,cellsCrch.size()); 
			tCrch = output.t; cellsCrchData=output.d; itCrch +=1;
			tnext = min(tnext,tCrch);
			}
		for (size_t i=0; i<cellsHrch.size(); i++) 
			{ 
			i1=clayRch[i];v=mesh.V()[i1];
			if ((CSu[i1]==0)&&(qHrch[i]>0)) { //recharge occurs only where no ghb, riv or drains + rch<0=evapo and no solutes evporaiton
				CSu[i1] = -qHrch[i]/v*cellsCrchData[i*4+2];
				}
			//Info<<"rch "<< i<<" nc " <<clayRch[i]<<" q "<<qHrch[i]<<" conc "<<cellsCrchData[i*4+2]<<endl;
			} 
	} //Info<<"Hrdh set"<<endl;
//Info<<" Rf "<<Rf<<" thetaw "<<thetaw<<endl;	
// ----------------------  equation  --------------------------
C.storePrevIter();//Info<<"c 10740 avt "<<C()[10740]<<endl;
fvScalarMatrix CEqn
	(
		Rf * thetaw * fvm::ddt(C) + Rf * C * fvc::ddt(thetaw) 
		+ fvm::div(phiw, C, "div(phiw,C)") 
		- fvm::laplacian(Deffw, C, "laplacian(Deffw,C)")
		+ thetaw * fvm::Sp(lbdaw,C)
		== - CSu - fvm::SuSp(CSp, C)
	); 
//for cells of fixed conc (it is supposed that we don't put cfix at otuflow)
if (cellsCfix.size()>0)
	{
	if (time>tCfix) {
		output = getCbuffer("cfix",itCfix,cellsCfix.size()); 
		tCfix = output.t; cellsCfixData=output.d; itCfix +=1;
		tnext = min(tnext,tCfix);
		}
	for (i=0; i<cellsCfix.size(); i++) {CEqn.setReference(cellsCfix[i],cellsCfixData[i*4+2]);}
	}
	
CEqn.solve();

//calculating dC, removing negative C and compensating mass
a=a* 0;dC=0;mxC=0;double dff=0;int imax=0;
for (i=0;i<ncell;i++) {mxC=max(mxC,C()[i]);sC+=C()[i];}
C=max(C,a); // removing negative
for (i=0;i<ncell;i++) {
	dff = mag(C()[i]-C.prevIter()[i]);
	if (dff>dC) {dC=dff;imax=i;}
	sC1+=C()[i];
	} 
	Info<<"cell "<<i1<<"dC "<<dC<<"mxC "<<mxC<<endl;
C *= sC/(sC1+SMALL); // removing negative value increase artifially the total
dC = dC/(mxC+SMALL);Info<<"dC "<<dC<<endl;

//Info<<" dC "<<dC<<endl;
dtForC = dCmax/(max(dC,0)+SMALL)*runTime.deltaTValue(); 
//if (activateThermal==1) {dtForC = min(dtForC,dTmax/(max(dT1,0)+SMALL)*runTime.deltaTValue());}

Info<< "dt "<<runTime.deltaTValue()<<" dC "<<dC<<" dtForC " << dtForC << endl; 
newDeltaT = min(dtForC, newDeltaT);
