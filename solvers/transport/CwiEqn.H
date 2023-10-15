dimensionedScalar USMALL("",dimLength/dimTime,1e-36);
// if (flowType==2) {thetaw = sw;} else {thetaw = eps*sw;}  // marche pas
thetaw = eps*sw;
//Deffw = (alphaT * mag(Ut))*tensor::I + (alphaL-alphaT) * (Ut*Ut)/(mag(Ut)+USMALL);
Deffw = (alphaT * mag(Uw))*tensor::I + (alphaL-alphaT) * (Uw*Uw)/(mag(Uw)+USMALL) + Dw0 * pow(thetaw,3.33)/pow(eps,2.)*tensor::I;
//Deffw = Dw0 * pow(thetaw,3.33)/pow(eps,2.)*tensor::I;
//forAll(Deffw(), i)  { if (sw[i]<sw_min[i]*1.5) {Deffw[i] *=1e2;} }; 

double conc,dC0,mxC,sC,sC1;int ic;size_t i;

float time = mesh.time().value()/day; //const Time&

//get the solution well data
if (cellsSwel.size()>0)
	{
	if (time>=tSwel) {
		output = getCbuffer("swel",itSwel,cellsSwel.size()); 
		tSwel = output.t; cellsSwelData=output.d; itSwel +=1;
		}
	for (i=0; i<cellsSwel.size(); i++) {condSwel[i] = static_cast<int>(cellsSwelData[i*4+2]);}
	}
//no solutions for drains
//get the solution ghb data
if (cellsSghb.size()>0)
	{
	if (time>=tSghb) {
		output = getCbuffer("sghb",itSghb,cellsSghb.size()); 
		tSghb = output.t; cellsSghbData=output.d; itSghb +=1;
		}
	for (i=0; i<cellsSghb.size(); i++) {condSghb[i] = static_cast<int>(cellsSghbData[i*4+2]);}
	}
// get fixed solutions 
if (cellsSfix.size()>0)
	{
	if (time>=tSfix) {
		output = getCbuffer("sfix",itSfix,cellsSfix.size()); 
		tSfix = output.t; cellsSfixData=output.d; itSfix +=1;
		}
	for (i=0; i<cellsSfix.size(); i++) {condSfix[i] = static_cast<int>(cellsSfixData[i*4+2]);}
	}
//get the solution recharge data
if (cellsSrch.size()>0)
	{
	if (time>=tSrch) {
		output = getCbuffer("srch",itSrch,cellsSrch.size()); 
		tSrch = output.t; cellsSrchData=output.d; itSrch +=1;
		}
	for (i=0; i<cellsSrch.size(); i++) 
		{
		condSrch[i] = static_cast<int>(cellsSrchData[i*4+2]);
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
		if (cellsHwel.size()>0)  // solution in well (even if no Swl data we ned to pump the conc)
		{
			for (i=0; i<cellsHwel.size(); i++) 
				{ 
				i1=cellsHwel[i];v=mesh.V()[i1];
				if (qHwel[i]<0) {CiSp[i1] = -qHwel[i]/v;} //in Sp term as discharge prop to local concentration
				else { CiSu[i1] = -qHwel[i]/v*solu_conc[condSwel[i]*ph_ncomp+ic]; } 
				}
		} //Info<<"Cweli set"<<endl;
		// solu at fixed head places (Sfix) need to be the same zones as Hfix (as we use Hfix[i])
		//Sfix are done with set reference below
		/*
		if (cellsSfix.size()>0)
		{
			for (i=0; i<cellsHfix.size(); i++) 
				{ //Info<<"i "<<i<<" condSfx "<<condSfx[i]<<endl;
				i1=cellsHfix[i];v=mesh.V()[i1];
				if (qHfix[i]<0) {CiSp[i1] = -qHfix[i]/v;} //in Sp term as discharge prop to local concentration
				else { CiSu[i1] = -qHfix[i]/v*solu_conc[condSfix[i]*ph_ncomp+ic]; } 
				} 
		} //Info<<"Cfxi set"<<endl;*/
		// conc at ghb places
		if (cellsHghb.size()>0)
		{
			for (i=0; i<cellsHghb.size(); i++) 
				{ 
				i1=cellsHghb[i];v=mesh.V()[i1];
				if (qHghb[i]<0) {CiSp[i1] = -(qHghb[i])/v;} //in Sp term as discharge prop to local concentration
				else {CiSu[i1] = -qHghb[i]/v*solu_conc[condSghb[i]*ph_ncomp+ic];}
				} 
		} //Info<<"Cghbi set"<<endl;
		//remove conc at drn places
		if (cellsHdrn.size()>0)
		{
			for (i=0; i<cellsHdrn.size(); i++) 
				{ i1=cellsHdrn[i];v=mesh.V()[i1];CiSp[i1] = -qHdrn[i]/v; } 
		} //Info<<"Cdrni set"<<endl;
		//recharge, for recharge however there shall be a solution at all surface cells so Srch and Hrch correspond
		if (cellsSrch.size()>0)
		{
			for (i=0; i<cellsSrch.size(); i++) 
				{ 
				v=mesh.V()[clayRch[i]];
				if (CiSu[clayRch[i]]==0) {
					CiSu[clayRch[i]] += -qHrch[i]/v*solu_conc[condSrch[i]*ph_ncomp+ic]; 
					} // recharge does not occur at drains, ghb, fix...
				//if (ic==0) {Info<<"rch "<< i<<" ilay " <<clayRch[i]<<" q "<<qHrch[i]<<" conc "<<solu_conc[condSrch[i]*ph_ncomp+ic]<<endl;}
				} 
		} //Info<<"Crchi set"<<endl;

		// for fixed solution cells with a search of the presence of fixed head or ghb cells
		icount = 0;float ccnd0=1e-6;float ccnd =0;int iH=1e5;int iG=1e5;
		//int ic0=2657;
		
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

		for (i=0; i<cellsSfix.size(); i++) {CwiEqn.setReference(cellsSfix[i],solu_conc[condSfix[i]*ph_ncomp+ic]);}
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
	} // of species loop
	
