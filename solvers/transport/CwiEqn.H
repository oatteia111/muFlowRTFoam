dimensionedScalar USMALL("",dimLength/dimTime,1e-36);
// if (flowType==2) {thetaw = sw;} else {thetaw = eps*sw;}  // marche pas
thetaw = eps*sw;
//Deffw = (alphaT * mag(Ut))*tensor::I + (alphaL-alphaT) * (Ut*Ut)/(mag(Ut)+USMALL);
Deffw = (alphaT * mag(Uw))*tensor::I + (alphaL-alphaT) * (Uw*Uw)/(mag(Uw)+USMALL) + Dw0 * pow(thetaw,3.33)/pow(eps,2.)*tensor::I;
//Deffw = Dw0 * pow(thetaw,3.33)/pow(eps,2.)*tensor::I;
//forAll(Deffw(), i)  { if (sw[i]<sw_min[i]*1.5) {Deffw[i] *=1e2;} }; 

double conc;int ic;size_t i;
mxCv.resize(ph_ncomp-4,0.);

time = mesh.time().value()/day; //nb tnext is the future time that a new condition will appear

//get the solution well data
if (cellsSwel.size()>0)
	{
	if (time>=tSwel) {
		output = getCbuffer("swel",itSwel,cellsSwel.size()); 
		tSwel = output.t; cellsSwelData=output.d; itSwel +=1;
		}
	for (i=0; i<cellsSwel.size(); i++) {condSwel[i] = static_cast<int>(cellsSwelData[i*4+2]);}
	//Info<<" tSwell "<<tSwel<<" cond 0 "<<condSwel[0]<<endl;
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
// get fixed solutions at h fix places
if (cellsShfix.size()>0)
	{
	if (time>=tShfix) {
		output = getCbuffer("shfix",itShfix,cellsShfix.size()); 
		tShfix = output.t; cellsShfixData=output.d; itShfix +=1;
		}
	for (i=0; i<cellsShfix.size(); i++) {condShfix[i] = static_cast<int>(cellsShfixData[i*4+2]);}
	//Info<<" tShfix "<<tShfix<<" cond 0 "<<condShfix[0]<<endl;
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
	for (i=0; i<cellsSrch.size(); i++) {condSrch[i] = static_cast<int>(cellsSrchData[i*4+2]);}
	}

float dC=0;dtForC=1e12;
/*
		for (j=0; j<nxyz;j++) {
			float c0=Cw[6]()[j];//if (j<2) {std::cout << "c0 1 "<<c0<<" ";}
			if ((c0>4e-4)&&(c0<6e-4)) {std::cout << "c "<<j<<" ";for (i=0; i<ph_ncomp;i++){std::cout<< Cw[i]()[j]<<" ";} std::cout<<"\n";} 
			sumT+=Cw[5]()[j];sumPF += c0;sumPFa +=Cw[7]()[j];
			}
		Info<< " sumPF  "<<sumT<<" "<<sumPF<<" "<<sumPFa << endl; 
		*/
/*
std::vector<double>sumC0(ph_ncomp-4,0.);
for (ic=4;ic<ph_ncomp;ic++) {for (j=0;j<nxyz;j++) {sumC0[ic-4]+=Cw[ic]()[j];} }
std::vector<double>sumC(ph_ncomp-4,0.);
*/
// main loop over species
unity=unity*0+1;
for (j=0;j<nxyz;j++) {Cw[0]()[j]=solu_conc[0];}  // water is alwas set to 55.34 from solution 0
for (ic=1; ic<ph_ncomp;ic++)
    {
		//Info<<"immob ic "<<immobile[ic]<<endl;
		volScalarField& Cwi = Cw[ic]();
		word nameCi ("Cw" + name(ic));
		if (immobile[ic]==1) {continue;}
		//Cw[ic]().storePrevIter();
		//Info<<"start Cwi "<<endl;
		//icnt=0;
		CiSu*=0; CiSp*=0; //previous cells with values must be destroyed
		// --------------------- internal conditions ------------------------
		if (cellsHwel.size()>0)  // solution in well (even if no Swl data we ned to pump the conc)
		{
			for (i=0; i<cellsHwel.size(); i++) 
				{ 
				i1=cellsHwel[i];v=mesh.V()[i1];unity[i1]=0;
				if (qHwel[i]<0) {CiSp[i1] = -qHwel[i]*(1+1e-6)/v;} //in Sp term as discharge prop to local concentration
				else { CiSu[i1] = -qHwel[i]/v*solu_conc[condSwel[i]*ph_ncomp+ic]; } 
				Info<< ic << " Csu wel "<<solu_conc[condSwel[i]*ph_ncomp+ic]<<" "<<CiSu[i1]<<" "<<Cwi[i1]<<endl;
				}
		} //Info<<"Cweli set"<<endl;
		// conc at ghb places
		if (cellsHghb.size()>0)
		{
			for (i=0; i<cellsHghb.size(); i++) 
				{ 
				i1=cellsHghb[i];v=mesh.V()[i1];
				if (qHghb[i]<0) {CiSp[i1] = -qHghb[i]/v;} //in Sp term as discharge prop to local concentration
				else {CiSu[i1] = -qHghb[i]/v*solu_conc[condSghb[i]*ph_ncomp+ic];}
				} 
		} //Info<<"Sghbi set"<<endl;
		// for fixed head cells
		if (cellsHfix.size()>0)
		{
			for (i=0; i<cellsHfix.size(); i++) 
				{ 
				i1=cellsHfix[i];v=mesh.V()[i1];//
				if (qHfix[i]<0) {CiSp[i1] = -qHfix[i]/v;} // outflow discharge prop to local concentration
				else {CiSu[i1] = -qHfix[i]/v*solu_conc[condShfix[i]*ph_ncomp+ic];}
				Info<< ic << " Csu hfix "<<solu_conc[condShfix[i]*ph_ncomp+ic]<<" "<<CiSu[i1]<<endl;
				}
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
				v=mesh.V()[clayRch[i]];
				if (CiSu[clayRch[i]]==0) {
					CiSu[clayRch[i]] += -qHrch[i]/v*solu_conc[condSrch[i]*ph_ncomp+ic]; 
					} // recharge does not occur at drains, ghb, fix...
				//if (ic==0) {Info<<"rch "<< i<<" ilay " <<clayRch[i]<<" q "<<qHrch[i]<<" conc "<<solu_conc[condSrch[i]*ph_ncomp+ic]<<endl;}
				} 
		} //Info<<"Srchi set"<<endl;
		
		//Rf = plugCwi.setRf(mesh,ic,Cwi); // this is the plugin for sorption,needs to be here (compn specific)
		//Info<<"Rf "<<ic<<" "<<Rf<<endl;
		fvScalarMatrix CwiEqn
		(		  
            Rf * thetaw * fvm::ddt(Cwi) + Rf * Cwi * fvc::ddt(thetaw)
            + fvm::div(phiw, Cwi, "div(phiw,C)")
            - fvm::laplacian(Deffw, Cwi, "laplacian(Deffw,C)") //*unity ne pas mettre thetaw
			+ CiSu + fvm::SuSp(CiSp, Cwi)
		);

		for (i=0; i<cellsSfix.size(); i++) {CwiEqn.setReference(cellsSfix[i],solu_conc[condSfix[i]*ph_ncomp+ic]);}

		CwiEqn.solve();

		//calculating dC, and time step (only for classical species)
		if (ic>3) {
			dC=0;mxC=0;dff=0;
			for (j=0;j<nxyz;j++) {mxC=max(mxC,Cwi[j]);}
			for (j=0;j<nxyz;j++) {
				dff = mag(Cwi[j]-Cw[ic]().prevIter()[j]);//Info<<"dff "<<j<<" "<<dff<<endl;
				dC= max(dC,dff);
				//sumC[ic-4]+=Cwi[j];
				} 
			dC = dC/(mxC+SMALL);Info<<"dC "<<dC<<" dtval " << runTime.deltaTValue()<<" c "<<Cwi[cellsHwel[0]]<<endl;
			dtForC = min(dtForC,dCmax/(max(dC,0)+SMALL)*runTime.deltaTValue()); 
		}
	}

Info<< "dt "<<runTime.deltaTValue()<<" dC "<<dC<<" dtForC " << dtForC << endl; 
//for (ic=0;ic<ph_ncomp-4;ic++) {Info<<ic<<" "<<sumC0[ic]<<" "<<sumC[ic]<<" "<<sumC[ic]-sumC0[ic]<<endl;}
newDeltaT = max(1e-4,min(dtForC, 1.2*runTime.deltaTValue()));
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
