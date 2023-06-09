dimensionedScalar USMALL("",dimLength/dimTime,1e-36);
// if (flowType==2) {thetaw = sw;} else {thetaw = eps*sw;}  // marche pas
thetaw = eps*sw;
//Deffw = (alphaT * mag(Ut))*tensor::I + (alphaL-alphaT) * (Ut*Ut)/(mag(Ut)+USMALL);
Deffw = (alphaT * mag(Uw))*tensor::I + (alphaL-alphaT) * (Uw*Uw)/(mag(Uw)+USMALL) + Dw0 * pow(thetaw,3.33)/pow(eps,2.)*tensor::I;
//Deffw = Dw0 * pow(thetaw,3.33)/pow(eps,2.)*tensor::I;
//forAll(Deffw(), i)  { if (sw[i]<sw_min[i]*1.5) {Deffw[i] *=1e2;} }; 

double conc,dC0,mxC,sC,sC1;int ic;

float time = mesh.time().value(); //const Time&

//get the solution well data
icount=0;
if (cellsSwl.size()>0)
	{
	float tnow = cellsHwData[2]*day; //time in files is in days
	while (time>=tnow) {icount += cellsSwl.size(); tnow = cellsSwlData[2+icount*4]*day; }
	if (icount>0) {icount -= cellsSwl.size();};
	for (int i=0; i<cellsSwl.size(); i++) {condSwl[i] = (int)cellsSwlData[2+(icount+i)*4+2];}
	}
//no solutions for drains
//get the solution ghb data
icount=0;
if (cellsSghb.size()>0)
	{
	float tnow = cellsSghbData[2]*day; //time in files is in days
	while (time>=tnow) {icount += cellsSghb.size(); tnow = cellsSghbData[2+icount*4]*day; }
	if (icount>0) {icount -= cellsSghb.size();};
	for (int i=0; i<cellsSghb.size(); i++) {condSghb[i] = (int)cellsSghbData[2+(icount+i)*4+2];}
	}
// get fixed solutions 
icount=0;
if (cellsSfx.size()>0)
	{
	float tnow = cellsSfxData[2]*day; //time in files is in days
	while (time>=tnow) {icount += cellsSfx.size(); tnow = cellsSfxData[2+icount*4]*day; }
	if (icount>0) {icount -= cellsSfx.size();};
	for (int i=0; i<cellsSfx.size(); i++) {condSfx[i] = (int)cellsSfxData[2+(icount+i)*4+2];}
	}
//get the solution recharge data
icount=0;
if (cellsSrch.size()>0)
	{
	float tnow = cellsSrchData[2]*day; 
	while (time>=tnow) {icount += cellsSrch.size(); tnow = cellsSrchData[2+icount*4]*day; }
	if (icount>0) {icount -= cellsSrch.size();};
	for (int i=0; i<cellsSrch.size(); i++) 
		{
		i1 = cellsSrch[i];
		for (ilay=0;ilay<nlay;ilay++) { if (sw[i1-ilay*ncell_lay]>0.999) {break;} ;} //search for the first confined layer from the top
		ilay1 = max(0, ilay-1); // ilay1 from top
		condSrch[i] = (int)cellsSrchData[2+(icount+i)*4+2];
		}
	}

//for (ilay=0;ilay<nlay;ilay++) {ic=ilay*ncell_lay+823;Info<<"cell 823 "<<ic<<" cwat "<<Cw[0]()[ic]<<endl;}

for (j=0;j<nxyz;j++) {Cw[0]()[j]=solu_conc[0];}  // water is alwas set to 55.34 from solution 0
dC=0;
for (ic=1; ic<ph_ncomp;ic++) //!!!!!!!!! no transport of water
    {
		volScalarField& Cwi = Cw[ic]();
		word nameCi ("Cw" + name(ic));
		Cwi.storePrevIter();
		//Info<<"start Cwi "<<endl;
		//set the solution to background for dry cells
		//forAll(Cwi(), celli) { if (sw[celli]<sw_min[i]*1.5) {Cwi[celli] = solu_conc[ic];} }
		icnt=0;
		CiSu*=0; CiSp*=0; //previous cells with values must be destroyed
		// --------------------- internal conditions ------------------------
		if (cellsHw.size()>0)  // solution in well (even if no Swl data we ned to pump the conc)
		{
			for (int i=0; i<cellsHw.size(); i++) 
				{ 
				i1=cellsHw[i];v=mesh.V()[i1];
				if (qHw[i]<0) {CiSp[i1] = -qHw[i]/v;} //in Sp term as discharge prop to local concentration
				else { CiSu[i1] = -qHw[i]/v*solu_conc[condSwl[i]*ph_ncomp+ic]; } 
				}
		} //Info<<"Hwi set"<<endl;
		// conc at fixed head places 
		if (cellsHfx.size()>0)
		{
			for (int i=0; i<cellsHfx.size(); i++) 
				{ //Info<<"i "<<i<<" condSfx "<<condSfx[i]<<endl;
				i1=cellsHfx[i];v=mesh.V()[i1];
				if (qHfx[i]<0) {CiSp[i1] = -qHfx[i]/v;} //in Sp term as discharge prop to local concentration
				else { CiSu[i1] = -qHfx[i]/v*solu_conc[condSfx[i]*ph_ncomp+ic]; } 
				} 
		} //Info<<"Hfx set"<<endl;
		// conc at ghb places
		if (cellsHghb.size()>0)
		{
			for (int i=0; i<cellsHghb.size(); i++) 
				{ 
				i1=cellsHghb[i];v=mesh.V()[i1];
				if (qHghb[i]<0) {CiSp[i1] = -(qHghb[i])/v;} //in Sp term as discharge prop to local concentration
				else {CiSu[i1] = -qHghb[i]/v*solu_conc[condSghb[i]*ph_ncomp+ic];}
				} 
		} //Info<<"Hghb set"<<endl;
		//remove conc at drn places
		if (cellsHdrn.size()>0)
		{
			for (int i=0; i<cellsHdrn.size(); i++) 
				{ i1=cellsHdrn[i];v=mesh.V()[i1];CiSp[i1] = -qHdrn[i]/v; } 
		} //Info<<"Hdrn set"<<endl;
		//recharge
		if (cellsSrch.size()>0)
		{
			for (int i=0; i<cellsSrch.size(); i++) 
				{ 
				v=mesh.V()[clayRch[i]];
				if (CiSu[clayRch[i]]==0) {
					CiSu[clayRch[i]] += -qHrch[i]/v*solu_conc[condSrch[i]*ph_ncomp+ic]; 
					} // recharge does not occur at drains, ghb, fix...
				//if (ic==0) {Info<<"rch "<< i<<" ilay " <<clayRch[i]<<" q "<<qHrch[i]<<" conc "<<solu_conc[condSrch[i]*ph_ncomp+ic]<<endl;}
				} 
		} //Info<<"Hrdh set"<<endl;

		// for fixed solution cells with a search of the presence of fixed head or ghb cells
		icount = 0;float ccnd0=1e-6;float ccnd =0;int iH=1e5;int iG=1e5;
int ic0=2657;
		
		/*
		if (cellsSfx.size()>0)
		{
			for (int i=0; i<cellsSfx.size(); i++) 
				{ 
				conc = solu_conc[condSfx[i]*ph_ncomp+ic]; //if ((i<10)&&(ic==0)) {Info<<"i "<<i<<" cSfx "<<condSfx[i]<<" conc "<<conc<<endl;}//conc = cellsCfxData[2+(icount+i)*4+2];
				i1=cellsSfx[i];v=mesh.V()[i1];ccnd = ccnd0*v/thk[i1];
				auto iter = std::find(cellsHfx.begin(), cellsHfx.end(), i1); // search in hFx cells
				iH = {std::distance(cellsHfx.begin(), iter)};
				auto iter1 = std::find(cellsHghb.begin(), cellsHghb.end(), i1); // search in ghb cells
				iG = {std::distance(cellsHghb.begin(), iter1)};
				if ((iH<cellsHfx.size()) && (qHfx[iH]>=0)) {CiSu[i1] = -qHfx[iH]/v*conc; }
				else if ((iG<cellsHghb.size()) && (qHghb[iG]>=0)) {CiSu[i1] = -qHghb[iG]/v*conc; }
				else if ((iH>=cellsHfx.size()) && (iG>=cellsHghb.size())) {CiSu[i1] = -ccnd*conc;CiSp[i1]=ccnd; }
				//Info<<"cfxG "<<i1<<" iG "<<iG<<" q "<<qHghb[iG]<<" conc "<<conc<<" Su "<<CiSu[i1]<<" Sp " <<CiSp[i1]<<endl;
				//if (i1 % ncell_lay == ic0) {Info<<"cfxH "<<i1<<" iH "<<iH<<" q "<<qHfx[iH]<<" conc "<<conc<<" Su "<<CiSu[i1]<<" Sp "<<CiSp[i1]<<endl;}
				//Info<<"cfx "<<i1<<" conc "<<conc<<" Su "<<CiSu[i1]<<" Sp "<<CiSp[i1]<<endl;}
				}
		}*/
		
		fvScalarMatrix CwiEqn
		(		  
            thetaw * fvm::ddt(Cwi) + Cwi * fvc::ddt(thetaw)
            + fvm::div(phiw, Cwi, "div(phiw,C)")
            - fvm::laplacian(Deffw, Cwi, "laplacian(Deffw,C)") // ne pa smettre thetaw
			+ CiSu + fvm::SuSp(CiSp, Cwi)
		);

		for (int i=0; i<cellsSfx.size(); i++) {CwiEqn.setReference(cellsSfx[i],solu_conc[condSfx[i]*ph_ncomp+ic]);}
		/*int nc;
		for (i=0;i<ncell_lay;i++) //set dry cells to solution 0
			{
			for (ilay=1;ilay<nlay;ilay++) 
				{ 
				nc = i+ilay*ncell_lay;//if (i==823) {Info<<"set ref 823 "<<nc<<" up "<<h()[nc]<<" down "<<h()[nc-ncell_lay]<<endl;}
				if (sw[nc]<sw_min[nc]*1.1) 
					{CwiEqn.setReference(nc,solu_conc[ic]);} 
				}
			}
		*/

		CwiEqn.solve();

		//dC calculation (not sure that it is good)
		/*if (ic>3) //not used for water, O, H,
		{
			dC0=0;sC=0;sC1=0;mxC=0;
			for (j=0;j<ncell;j++) {mxC=max(mxC,Cwi()[j]);}
			//Cwi=max(Cwi,Cwi*0); // removing negative
			for (j=0;j<ncell;j++) {dC0= max(dC0,Cwi()[j]-Cwi.prevIter()[j]);}
			dC0 = dC0/(mxC+SMALL); // dC is relative
			Info <<"comp "<<ic<<" dC0 "<<dC0<<" dC "<<dC<<endl;
		}
		dC = max(dC,dC0); // we consider the species that leads to the max dC*/
	}


/* int nc;float minCwat=100.;
for (i=0;i<ncell_lay;i++) //set cw_lay i = cw_lay i-1 if sw=sw_min
	{
	for (ilay=1;ilay<nlay;ilay++) 
		{ 
		nc = i+ilay*ncell_lay;if (Cw[0]()[nc]<minCwat) {imin=nc;minCwat=Cw[0]()[nc];}
		}
	}
Info<<"minCw "<<imin<<" "<<minCwat<<" in frk "<<freak.c[imin]<<endl; */
