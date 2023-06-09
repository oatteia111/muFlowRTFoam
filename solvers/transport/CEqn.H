
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

// --------------------- internal conditions ------------------------
	float time = mesh.time().value(); //const Time&
	CSu*=0; CSp*=0; a1=1; double conc;//previous cells with values must be destroyed
				
	icount = 0;
	if (cellsHw.size()>0)
	{
		//get the concentration well data
		float tnow = cellsHwData[2]*day; //time in files is in days
		while (time>=tnow) {icount += cellsCw.size(); tnow = cellsCwData[2+icount*4]*day; }
		if (icount>0) {icount -= cellsCw.size();};
		for (int i=0; i<cellsCw.size(); i++) 
			{ 
			i1=cellsCw[i];v=mesh.V()[i1];
			if (qHw[i]<0) {CSp[i1] = -qHw[i]/v;} //in Sp term as discharge prop to local concentration, a positive Sp will make solution implicit
			else {CSu[i1] = -qHw[i]/v*cellsCwData[2+(icount+i)*4+2];} //the injected mass is the water flow x the conc
			//Info<<"cellw "<<cellsCw[i]<<" q "<<qHw[i]<<" cloc "<<cloc<<" Su "<<CSu[cellsCw[i]]<<" Sp "<<CSp[cellsCw[i]]<<endl;
			}
	}
	// for fixed head cells only for the outflow of conc (1st step) -> Hfx because Cghb may not be set here
	icount = 0;
	if (cellsHfx.size()>0)
	{
		float tnow = cellsHfxData[2]*day; //time in files is in days
		while (time>=tnow) {icount += cellsHfx.size(); tnow = cellsHfxData[2+icount*4]*day; }
		if (icount>0) {icount -= cellsHfx.size();};
		for (int i=0; i<cellsHfx.size(); i++) 
			{ 
			i1=cellsHfx[i];v=mesh.V()[i1];//
			if (qHfx[i]<0) {CSp[i1] = -qHfx[i]/v;} // outflow qHfx<0 discharge prop to local concentration
			else {CSu[i1] = -qHfx[i]/v*cellsCfxData[2+(icount+i)*4+2];}
			}
	}
	// for ghb cells only for the outflow of conc (1st step) -> Hghb because Cghb may not be set here
	icount = 0;
	if (cellsHghb.size()>0)
	{
		float tnow = cellsHghbData[2]*day; //time in files is in days
		while (time>=tnow) {icount += cellsHghb.size(); tnow = cellsHghbData[2+icount*4]*day; }
		if (icount>0) {icount -= cellsHghb.size();};
		for (int i=0; i<cellsHghb.size(); i++) 
			{ 
			i1=cellsHghb[i];v=mesh.V()[i1];//
			if (qHghb[i]<0) {CSp[i1] = -qHghb[i]/v;} // outflow discharge prop to local concentration
			else {CSu[i1] = -qHghb[i]/v*cellsCghbData[2+(icount+i)*4+2];}
			//Info<<"ghb "<< i<<" i1 " <<i1<<" q "<<qHghb[i]<<" conc "<<cellsCghbData[2+(icount+i)*4+2]<<endl;
			}
	}
	// for drain cells (only out), same as above, except that drains are constant in time
	icount = 0;
	if (cellsHdrn.size()>0)
	{
		for (int i=0; i<cellsHdrn.size(); i++) 
			{ 
			i1=cellsHdrn[i];v=mesh.V()[i1];//
			if (qHdrn[i]<0) {CSp[i1] = -qHdrn[i]/v;} // outflow qHdrn<0 discharge prop to local concentration
			}
	}
	// for cfx conc (if no hfix)
	icount = 0;
	if (cellsCfx.size()>0)
	{
		float tnow = cellsCfxData[2]*day; //time in files is in days
		while (time>=tnow) {icount += cellsCfx.size(); tnow = cellsCfxData[2+icount*4]*day; }
		if (icount>0) {icount -= cellsCfx.size();};
		for (int i=0; i<cellsCfx.size(); i++) {i1=cellsCfx[i]; condCfx[i] = cellsCfxData[2+(icount+i)*4+2];}
	}

	//recharge
	icount=0;int nc;
	if (cellsCrch.size()>0) // recharge is different, required only if C>0
	{
		float tnow = cellsCrchData[2]*day; //time in files is in days
		while (time>=tnow) {icount += cellsCrch.size(); tnow = cellsCrchData[2+icount*4]*day; }
		if (icount>0) {icount -= cellsCrch.size();};
		for (int i=0; i<cellsCrch.size(); i++) 
			{ 
			v=mesh.V()[clayRch[i]];
			if (CSu[clayRch[i]]==0) { //recharge occurs only where no ghb, riv or drains
				CSu[clayRch[i]] += -qHrch[i]/v*cellsCrchData[2+(icount+i)*4+2];
				}
			//else {Info<<"rch "<< i<<" nc " <<clayRch[i]<<" q "<<qHrch[i]<<" conc "<<cellsCrchData[2+(icount+i)*4+2]<<endl;}
			} 
	} //Info<<"Hrdh set"<<endl;


//Info<<"0 "<<i<<" csu "<<CSu[0]<<" csp "<<CSp[0]<<" c "<<C()[0]<<endl;
	
//get the fixed conc conditions
/*	if (cellsCfx.size()>0)
	{
		float tnow = cellsCfxData[2]*day; //time in files is in days
		while (time>=tnow) {icount += cellsCfx.size(); tnow = cellsCfxData[2+icount*4]*day; }
		if (icount>0) {icount -= cellsCfx.size();};
		for (int i=0; i<cellsCfx.size(); i++) {condCfx[i] = cellsCfxData[2+(icount+i)*4+2];} //Info<<"cfx "<<b[i]<<endl;} // discharge must be in s-1
	} */
//forAll(C(), i)  { if (sw[i]<sw_min[i]*1.5) {cellSet -noFlip -cellNumber $i -region interior;} } chatGPT -> wrong
//forAll(C(), i)  { if (sw[i]<sw_min[i]*1.5) {a1[i]=0;} }
	
// ----------------------  equation  --------------------------
C.storePrevIter();//Info<<"c 10740 avt "<<C()[10740]<<endl;
fvScalarMatrix CEqn
	(
		R * thetaw * fvm::ddt(C) + R * C * fvc::ddt(thetaw) 
		+ fvm::div(phiw, C, "div(phiw,C)") 
		- fvm::laplacian(Deffw, C, "laplacian(Deffw,C)")
		+ thetaw * fvm::Sp(lbdaw,C)
		+ CSu + fvm::SuSp(CSp, C)
	); 

for (int i=0; i<cellsCfx.size(); i++) {CEqn.setReference(cellsCfx[i],condCfx[i]);}
CEqn.solve();
//Info<<"c 10740 apr "<<C()[10740]<<endl;
//Info <<"resid "<<CwEqn.residual()()[67171]<<endl;//.value()<<endl;
//Info << "Ceqn solve, alphaL "<<alphaL<<" theta 5 "<<thetaw[5]<<" Deffw 5"<<Deffw[5]<<endl;
//scalar dCdTmax = max(mag(fvc::ddt(Cw))).value();Info<<"dCdTmax "<<dCdTmax<<endl;
//scalar dC = dCdTmax*runTime.deltaTValue()*1;
//Info<<"max dCw " <<max(mag(Cw-Cw.prevIter()))<<endl;

//calculating dC, removing negative C and compensating mass
a *= 0;dC=0;mxC=0;double dff=0;int i1=0;
for (i=0;i<ncell;i++) {mxC=max(mxC,C()[i]);sC+=C()[i];}
C=max(C,a); // removing negative
for (i=0;i<ncell;i++) {
	dff = mag(C()[i]-C.prevIter()[i]);
	if (dff>dC) {i1=i;}
	dC= max(dC,dff);sC1+=C()[i];
	} 
	Info<<"cell "<<i1<<"dC "<<dC<<"mxC "<<mxC<<endl;
C *= sC/(sC1+SMALL); // removing negative value increase artifially the total
dC = dC/(mxC+SMALL);Info<<"dC "<<dC<<endl;


