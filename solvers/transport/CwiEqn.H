const scalar dCmax = runTime.controlDict().lookupOrDefault("dCmax", 0.02);

dimensionedScalar USMALL("",dimLength/dimTime,1e-36);
// if (flowType==2) {thetaw = sw;} else {thetaw = eps*sw;}  // marche pas
thetaw = eps*sw;
//Deffw = (alphaT * mag(Ut))*tensor::I + (alphaL-alphaT) * (Ut*Ut)/(mag(Ut)+USMALL);
Deffw = (alphaT * mag(Uw))*tensor::I + (alphaL-alphaT) * (Uw*Uw)/(mag(Uw)+USMALL) + Dw0 * pow(thetaw,3.33)/pow(eps,2.)*tensor::I;
//Deffw = Dw0 * pow(thetaw,3.33)/pow(eps,2.)*tensor::I;
//forAll(Deffw(), i)  { if (sw[i]<sw_min[i]*1.5) {Deffw[i].zz()=0;} }; 

double conc;int ic;size_t i;int jmax;
mxCv.resize(ph_ncomp-4,0.);
bcCwi.resize(ncell,0);

time = mesh.time().value() ; //nb tnext is the future time that a new condition will appear

if (flowType<4) 
{
	//no solutions for drains
	//get the solution ghb data
	std::cout<<"nb cell ghb "<<cellsSghb.size()<<" cond "<<condSghb.size()<<" it "<< itSghb<<" data "<<cellsSghbData.size()<<" bc "<<bcCwi.size() <<"\n";
	if (cellsSghb.size()>0)
		{
		if (time>=tSghb) {
			output = getCbuffer("sghb",itSghb,cellsSghb.size()); 
			tSghb = output.t; cellsSghbData=output.d; itSghb +=1;
			tnext = min(tnext,tSghb);
			}
		for (i=0; i<cellsSghb.size(); i++) 
			{condSghb[i] = static_cast<int>(cellsSghbData[i*4+2]);i1=cellsSghb[i];bcCwi[i1]=1;}
		Info<<" tSghb "<<tSghb<<" cond 0 "<<condSghb[0]<<endl;
		}
	// get fixed solutions at h fix places
	if (cellsShfix.size()>0)
		{
		if (time>=tShfix) {
			output = getCbuffer("shfix",itShfix,cellsShfix.size()); 
			tShfix = output.t; cellsShfixData=output.d; itShfix +=1;
			tnext = min(tnext,tShfix);
			}
		for (i=0; i<cellsShfix.size(); i++) 
			{condShfix[i] = static_cast<int>(cellsShfixData[i*4+2]);i1=cellsShfix[i];bcCwi[i1]=1;}
		Info<<" tShfix "<<tShfix<<" cond 0 "<<condShfix[0]<<endl;
		}
	//get the solution recharge data
	if (cellsSrch.size()>0)
		{
		if (time>=tSrch) {
			output = getCbuffer("srch",itSrch,cellsSrch.size()); 
			tSrch = output.t; cellsSrchData=output.d; itSrch +=1;
			tnext = min(tnext,tSrch);
			}
		for (i=0; i<cellsSrch.size(); i++) 
			{condSrch[i] = static_cast<int>(cellsSrchData[i*4+2]);i1=cellsSrch[i];bcCwi[i1]=1;}
		Info<<" tSrch "<<tSrch<<" cond 0 "<<condSrch[0]<<endl;
		}
} // end flowtype<4
else 
{
	if (cellsSpfix.size()>0)
	{
	if (time>=tSpfix) {
		output = getCbuffer("spfix",itSpfix,cellsSpfix.size()); 
		tSpfix = output.t; cellsSpfixData=output.d; itSpfix +=1;
		tnext = min(tnext,tSpfix);
		}
	for (i=0; i<cellsSpfix.size(); i++) 
		{condSpfix[i] = static_cast<int>(cellsSpfixData[i*4+2]);i1=cellsSpfix[i];bcCwi[i1]=1;}
	Info<<" tSpfix "<<tSpfix<<" cond 0 "<<condSpfix[0]<<endl;
	}
}

// get fixed solutions 
if (cellsSfix.size()>0)
	{
	if (time>=tSfix) {
		output = getCbuffer("sfix",itSfix,cellsSfix.size()); 
		tSfix = output.t; cellsSfixData=output.d; itSfix +=1;
		tnext = min(tnext,tSfix);
		}
	for (i=0; i<cellsSfix.size(); i++) 
		{condSfix[i] = static_cast<int>(cellsSfixData[i*4+2]);i1=cellsSfix[i];bcCwi[i1]=1;}
	}
//get the solution well data
if (cellsSwel.size()>0)
	{
	if (time>=tSwel) {
		output = getCbuffer("swel",itSwel,cellsSwel.size()); 
		tSwel = output.t; cellsSwelData=output.d; itSwel +=1;
		tnext = min(tnext,tSwel);
		}
	for (i=0; i<cellsSwel.size(); i++) 
		{condSwel[i] = static_cast<int>(cellsSwelData[i*4+2]);i1=cellsSwel[i];bcCwi[i1]=1;}
	Info<<" tSwell "<<tSwel<<" cond 0 "<<condSwel[0]<<endl;
	}

float dC=0;dtForC=1e12;

// main loop over species
unity=unity*0+1;

for (j=0;j<nxyz;j++) {Cw[0]()[j]=solu_conc[0];}  // water is alwas set to 55.34 from solution 0
for (ic=1; ic<ph_ncomp;ic++)
    {
		//Info<<"immob ic "<<immobile[ic]<<endl;
		volScalarField& Cwi = Cw[ic]();
		word nameCi ("Cw" + name(ic));
		if (immobile[ic]==1) {continue;}
		//for (j=0;j<ncell-nxyz;j++) {Cw[ic]()[rinactive[j]] = solu_conc[ic];} // set inacive cells to background
		//for (j=0;j<ncell;j++) {if (sw[j]<1e-3) {Cw[ic]()[j]=solu_conc[ic];}
			//if (j<10) {Info<<"dry cell "<<j<<" sw "<<sw[j]<<" "<<Cw[ic]()[j]<<" Csu "<<CiSu[j]<<endl;}
			//} //set dry cells to background
		//icnt=0;
		CiSu*=0; CiSp*=0; //previous cells with values must be destroyed
		// --------------------- internal conditions ------------------------
		if (flowType<4)
		{
			if (cellsHwel.size()>0)  // solution in well (even if no Swl data we ned to pump the conc)
			{
				for (i=0; i<cellsHwel.size(); i++) 
					{ 
					i1=cellsHwel[i];v=mesh.V()[i1];unity[i1]=0;
					if (qHwel[i]<0) {CiSp[i1] = -qHwel[i]*(1+1e-6)/v;} //in Sp term as discharge prop to local concentration
					else { CiSu[i1] = -qHwel[i]/v*solu_conc[condSwel[i]*ph_ncomp+ic]; } 
					}
				Info<< ic << " Csu wel "<<solu_conc[condSwel[0]*ph_ncomp+ic]<<" "<<CiSu[i1]<<" "<<Cwi[i1]<<endl;
			} //Info<<"Cweli set"<<endl;
			// conc at ghb places
			if (cellsHghb.size()>0)
			{
				for (i=0; i<cellsHghb.size(); i++) 
					{ 
					i1=cellsHghb[i];v=mesh.V()[i1];
					if (qHghb[i]<0) {CiSp[i1] = -qHghb[i]/v;} //in Sp term as discharge prop to local concentration
					else {CiSu[i1] = -qHghb[i]/v*solu_conc[condSghb[i]*ph_ncomp+ic];}
					//else {CiSu[i1] = -solu_conc[condSghb[i]*ph_ncomp+ic];CiSp[i1]=1;} //does not seem to work properly
					} 
				//Info<< ic << " Csu ghb "<<solu_conc[condSghb[0]*ph_ncomp+ic]<<" "<<CiSu[i1]<<" "<<Cwi[i1]<<endl;
			} //Info<<"Sghbi set"<<endl;
			// for fixed head cells
			if (cellsHfix.size()>0)
			{
				for (i=0; i<cellsHfix.size(); i++) 
					{ 
					i1=cellsHfix[i];v=mesh.V()[i1];//
					if (qHfix[i]<0) {CiSp[i1] = -qHfix[i]/v;} // outflow discharge prop to local concentration
					else {CiSu[i1] = -qHfix[i]/v*solu_conc[condShfix[i]*ph_ncomp+ic];}
					}
				//Info<< ic << " Csu hfix "<<solu_conc[condShfix[0]*ph_ncomp+ic]<<" "<<CiSu[i1]<<endl;
			}
			//remove conc at drn places
			if (cellsHdrn.size()>0)
			{
				for (i=0; i<cellsHdrn.size(); i++) 
					{ i1=cellsHdrn[i];v=mesh.V()[i1];CiSp[i1] = -qHdrn[i]/v; } 
			} //Info<<"Cdrni set"<<endl;
			//recharge, for recharge however there shall be a solution at all surface cells so Srch and Hrch correspond
			if (cellsHrch.size()>0)
			{
				for (i=0; i<cellsHrch.size(); i++) 
					{ 
					i1=clayRch[i];v=mesh.V()[i1];
					if ((CiSu[i1]==0)&&(qHrch[i]>0)) { // recharge does not occur at drains, ghb, fix; and if rch<0: evapo, no solutes
						CiSu[i1] += -qHrch[i]/v*solu_conc[condSrch[i]*ph_ncomp+ic]; 
						} 
					//if (i<5) {Info<<"rch "<< i<<" i1 " <<i1<<" q "<<qHrch[i]<<" conc "<<solu_conc[condSrch[i]*ph_ncomp+ic]<<endl;}
					} 
			} //Info<<"Srchi set"<<endl;
			// for dry cells, set conc as background
			if (flowType==2) 
				{for (j=0;j<ncell;j++) {if (sw[j]<sw_min[j]*1.1) {CiSu[j]=-solu_conc[ic];CiSp[j]=1;} } }

		} //end of flowtype<4 condition
		
		else // here flowType=4, conditions are for water
		{
		if (cellsPfix.size()>0)
			{
				for (i=0; i<cellsPfix.size(); i++) 
					{ 
					i1=cellsPfix[i];v=mesh.V()[i1];//Info<<"in cwi pfix "<<i<<" "<<i1<<" q "<<qPfix[i]<<endl;
					if (qPfix[i]<=0) {CiSp[i1] = -qPfix[i]/v;} // outflow discharge prop to local concentration
					else {CiSu[i1] = -qPfix[i]/v*solu_conc[condSpfix[i]*ph_ncomp+ic];}
					//if (i==0) {Info<< " Csu pfix (0) "<<solu_conc[condSpfix[i]*ph_ncomp+ic]<<" "<<CiSu[i1]<<endl;}
					}
			}
		if (cellsPwel.size()>0)  // solution in well (even if no Swl data we ned to pump the conc)
			{
				for (i=0; i<cellsPwel.size(); i++) 
					{ 
					i1=cellsPwel[i];v=mesh.V()[i1];unity[i1]=0;Info<<"qPwel "<<qPwel[i]<<"\n";
					if (qPwel[i]<=0) {CiSp[i1] = -qPwel[i]*(1+1e-6)/v;} //in Sp term as discharge prop to local concentration
					else {CiSu[i1] = -qFwel[i]/v*solu_conc[condSwel[i]*ph_ncomp+ic]; } 
					//if (i==0) {Info<<i<< " Csu wel "<<solu_conc[condSwel[i]*ph_ncomp+ic]<<" "<<CiSu[i1]<<" "<<Cwi[i1]<<endl;}
					}
			} //Info<<"Cweli set"<<endl;

		}
		//for (j=1;j<ncell-nxyz;j++) {CiSu[rinactive[j]]=0;CiSp[rinactive[j]]=0;}
		//Rf = plugCwi.setRf(mesh,ic,Cwi); // this is the plugin for sorption,needs to be here (compn specific)
		//Info<<"Rf "<<ic<<" "<<Rf<<endl;
		

		fvScalarMatrix CwiEqn
		(		  
            Rf * thetaw * fvm::ddt(Cwi) + Rf * Cwi * fvc::ddt(thetaw)
            + fvm::div(phiw, Cwi, "div(phiw,C)")
            - fvm::laplacian(Deffw, Cwi, "laplacian(Deffw,C)") //*unity ne pas mettre thetaw
			== - CiSu - fvm::SuSp(CiSp, Cwi)
		);

		for (i=0; i<cellsSfix.size(); i++) {CwiEqn.setReference(cellsSfix[i],solu_conc[condSfix[i]*ph_ncomp+ic]);}

		CwiEqn.solve();
		Info<<"Cw (0) "<<Cw[ic]()[0]<<endl;

		//calculating dC, and time step (only for classical species)
		if (ic>3) {
			dC=0;mxC=0;dff=0;jmax=0;
			for (j=0;j<nxyz;j++) {mxC=max(mxC,Cw[ic]()[ractive[j]]);} //Info<<"mxC "<<mxC<<endl;
			for (j=0;j<nxyz;j++) if (bcCwi[j]==0) {
				dff = mag(Cw[ic]()[ractive[j]]-Cw[ic]().prevIter()[ractive[j]]);
				if (dff>dC) {dC=dff;jmax=j;}
				} 
			dC = dC/(mxC+SMALL);
			dtForC = min(dtForC,dCmax/(max(dC,0)+SMALL)*runTime.deltaTValue()); 
			Info<<"dC "<<dC<<" mxC " << mxC<<" dtForC "<<dtForC<<" jmax "<<jmax<<" ract(j) "<<ractive[jmax]<<" Csu "<<CiSu[jmax]<<endl;
		}
	}

Info<< "dt "<<runTime.deltaTValue()<<" dC "<<dC<<" dtForC " << dtForC << endl; 
//for (ic=0;ic<ph_ncomp-4;ic++) {Info<<ic<<" "<<sumC0[ic]<<" "<<sumC[ic]<<" "<<sumC[ic]-sumC0[ic]<<endl;}
newDeltaT = max(minDeltaT,min(dtForC, newDeltaT));
		/*sumT=0;sumPF=0;sumPFa=0;
		for (j=0; j<nxyz;j++) {
			float c0=Cw[6]()[j];//if (j<2) {std::cout << "c0 1 "<<c0<<" ";}
			if ((c0>4e-4)&&(c0<6e-4)) {
				std::cout << "c "<<j<<" ";for (i=0; i<ph_ncomp;i++){std::cout<< Cw[i]()[j]<<" ";} std::cout<<"\n";
				} 
			sumT+=Cw[5]()[j];sumPF += c0;sumPFa +=Cw[7]()[j];
			}
		Info<< " sumPF  "<<sumT<<" "<<sumPF<<" "<<sumPFa << endl; 
	*/
