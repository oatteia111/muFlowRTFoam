dimensionedScalar USMALL("",dimLength/dimTime,1e-36);
// if (flowType==2) {thetaw = sw;} else {thetaw = eps*sw;}  // marche pas
thetaw = eps*sw;
//Deffw = (alphaT * mag(Ut))*tensor::I + (alphaL-alphaT) * (Ut*Ut)/(mag(Ut)+USMALL);
Deffw = (alphaT * mag(Uw))*tensor::I + (alphaL-alphaT) * (Uw*Uw)/(mag(Uw)+USMALL) + Dw0 * pow(thetaw,3.33)/pow(eps,2.)*tensor::I;
//Deffw = Dw0 * pow(thetaw,3.33)/pow(eps,2.)*tensor::I;
//forAll(Deffw(), i)  { if (sw[i]<sw_min[i]*1.5) {Deffw[i] *=1e2;} }; 

double conc;int ic;size_t i;
mxCv.resize(ph_ncomp-4,0.);

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

double dC=0;double sdC=0;

for (j=0;j<nxyz;j++) {Cw[0]()[j]=solu_conc[0];}  // water is alwas set to 55.34 from solution 0
for (ic=1; ic<ph_ncomp;ic++)
    {
		volScalarField& Cwi = Cw[ic]();
		word nameCi ("Cw" + name(ic));
		Cw[ic]().storePrevIter();
		//Info<<"start Cwi "<<endl;
		icnt=0;
		CiSu*=0; CiSp*=0; //previous cells with values must be destroyed
		// --------------------- internal conditions ------------------------
		if (cellsHwel.size()>0)  // solution in well (even if no Swl data we ned to pump the conc)
		{
			for (i=0; i<cellsHwel.size(); i++) 
				{ 
				i1=cellsHwel[i];v=mesh.V()[i1];
				if (qHwel[i]<0) {CiSp[i1] = -qHwel[i]*(1+1e-6)/v;} //in Sp term as discharge prop to local concentration
				else { CiSu[i1] = -qHwel[i]/v*solu_conc[condSwel[i]*ph_ncomp+ic]; } 
				}
		} //Info<<"Cweli set"<<endl;
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
		//icount = 0;float ccnd0=1e-6;float ccnd =0;int iH=1e5;int iG=1e5;
		
		fvScalarMatrix CwiEqn
		(		  
            thetaw * fvm::ddt(Cwi) + Cwi * fvc::ddt(thetaw)
            + fvm::div(phiw, Cwi, "div(phiw,C)")
            - fvm::laplacian(Deffw, Cwi, "laplacian(Deffw,C)") // ne pa smettre thetaw
			+ CiSu + fvm::SuSp(CiSp, Cwi)
		);

		for (i=0; i<cellsSfix.size(); i++) {CwiEqn.setReference(cellsSfix[i],solu_conc[condSfix[i]*ph_ncomp+ic]);}

		CwiEqn.solve();
	}

/*		//calculating dC, (only for true species)
		if (ic>3) 
			{
			double dff=0;int i1=0;double sC = 0;double dC0=0;double mxC0=0;
			//for (int i=0;i<ncell;i++) {mxC0=max(mxC0,Cwi.prevIter()[i]);} //strange, ifwe remove that dC0 is wrong?
			for (int i=0;i<ncell;i++) {mxCv[ic-4]=max(mxCv[ic-4],Cwi[i]);} //sC+=Cwi[i];
			if (mxCv[ic-4]>0) {//Cwi=max(Cwi,a); // removing negative
				for (int i=0;i<ncell;i++) {
					dff = mag(Cwi[i]-Cw[ic]().prevIter()[i])/(mxCv[ic-4]+SMALL);
					dC0= max(dC0,dff);//sC1+=Cwi[i]; //species max
					} 
				//Info<<"dC0 "<<dC0<<" mxC0 "<<mxC0<<" mxC "<<mxC<<endl;
				dC=max(dC0,dC);
				}
			}
		
	} // end of species loop
		
dtForC = dCmax/(max(dC,0)+SMALL)*runTime.deltaTValue(); 
//if (activateThermal==1) {dtForC = min(dtForC,dTmax/(max(dT1,0)+SMALL)*runTime.deltaTValue());}

Info<< "dt "<<runTime.deltaTValue()<<" dC "<<dC<<" dtForC " << dtForC << endl; 
newDeltaT = min(dtForC, 1.2*runTime.deltaTValue());
runTime.setDeltaT (min (newDeltaT,maxDeltaT) ); */