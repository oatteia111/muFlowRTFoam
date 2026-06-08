tcnt++;
if (rSteps<0) {if (tcnt>-rSteps-1) {tcnt=0;} }
if (rSteps>0) {if (rcnt>=rSteps) {rcnt=1;} }
int flg =0;
//if (activateReaction==1 && tcnt==rSteps-1)
if (rSteps>0) {if (runTime.value() >= wTimes[itwstep-1]+reactStep*rcnt) flg=1;} // here the time to write is a portion of current time period
if (rSteps<0) {if (tcnt == -rSteps-1) flg=1;} // here the time is a number of flow/transport time steps
Info<<"for react tcnt "<<tcnt<<" wtime "<<wTimes[itwstep-1]<<" lim "<<wTimes[itwstep-1]+reactStep*rcnt<<" rstep "<<rSteps<<" tcnt "<<tcnt<<" rcnt "<<rcnt<<" flg "<<flg<<endl;
if (activateReaction==1 && flg)
{
	rcnt++;flg=0;
	//if (tcnt>rSteps) {tcnt=0;}
	// ***************** find the cells where the conc has changed to calculate there
	//deltaTchem = transportProperties.lookupOrDefault<scalar>("deltaTchem",86400);Info<<"dtchem in reac "<<deltaTchem<<endl;
	std::vector<double> rchange(nxyz,1.);
	if (ncell>1000) {
		for (j=0; j<nxyz;j++) {
			rchange[j]=0.;
			int rj=ractive[j];
			for (i=4; i<ph_ncomp; i++)
			{
				//if (abs(Cw[i]().prevIter()[rj]-Cw[i]()[rj])/(Cw[i]()[rj]+1e-20)>1e-10) {rchange[j] = 1.;}
				if (Cw[i]()[rj]>1e-12) {rchange[j] = 1.;}
			} 
			for (i=0; i<ph_gcomp; i++)
			{
				if (abs(gm_ph[i*nxyz+j]-Cg[i]()[rj])/(gm_ph[i*nxyz+j]+1e-20)<1e-10 || sw[rj]<sw_min[rj]*1.01) {rchange[j] = 0.;} else {rchange[j] = 1.;}
			}
		}
	}
	
	// ********************** transfer to phreeqc for c and g, for g we send moles and not pressures, phreeqc does not know the cell volume it considers to be one
	//int icnt = 0;
	forAll(Cw,i) { for (j=0; j<nxyz;j++) {c_ph[i*nxyz+j] = Cw[i]()[ractive[j]];} } 
	forAll(Cg,i) {
			//Vmol[j] = 24.47*(273.15+T()[j])/(273.15+25.);
		for (j=0; j<nxyz;j++) { 
			gm_ph[i*nxyz+j] = max(Cg[i]()[j]*gvol[j],1e-16);// Cg in mol/Lgaz, gm in mol/RV
			} //Cg in moles/L and Vm in L/mol
		} 
	//for (j=0;j<6;j++) { for (i=4; i<ph_ncomp;i++){Info<<Cw[i]()[j]<<" ";} Info<<endl;}
	//Info<<" phqVm[0] "<<Vmol[0]<<" "<<endl;
	auto start = std::chrono::high_resolution_clock::now();
	
	//################# RUN PHREEQC   ################
	// set saturations using rchange
	t_ph.resize(nxyz); p_ph.resize(nxyz);// showul dnot be necessary but seems required???
	//Info<<"rchange ";
	for (j=0; j<nxyz;j++) {rchange[j] *= sw[j];t_ph[j]=T()[j]; } //if(j<10) {Info<<" "<<rchange[j];} } Info<<endl;
	//for (j=0; j<nxyz;j++) {p_ph[j]=p[j]/atmPa;}
	if (ph_gcomp>0) {freak.setGvol(gvol);} // set gas volume in phreeqc
	freak.setWsat(rchange); // rchange for the calculation domain, with 0 outside, sw saturation
	freak.setC(c_ph);//transfer c_ph to freak : it does not work to send directly to freak.c
	freak.setGm(gm_ph);//transfer gm_ph to freak
	if (activateThermal==1) {freak.setTemp(t_ph);}
	//freak.setP(p_ph);//transfer pressure to freak
	freak.setTstep((runTime.value()-oldTimeReac)*tunits); //Info<<" this tme "<< runTime.value()<<" old "<<presentTime<<endl;//the calculation time shall include all time since las phreeqc run
	Info << "running phreeqc dt "<<(runTime.value()-oldTimeReac)*tunits<<endl;
	int a0= phqRun(freak);//Info<<"end results "<<a0<<endl;
	//for (j=0;j<6;j++) { for (i=4; i<ph_ncomp;i++){Info<<freak.c[i*nxyz+j]<<" ";} Info<<endl;}
	if (a0==-7) { // irm _fail : what to do
		forAll(Cw,ic) { if (immobile[ic]==0) {Cw[ic]()=Cw[ic]().prevIter();} } // back to previous conc values
		tcnt=rSteps-2;//insure that phreeqc will be run on next time step !!!!! oldTimeReac to be corrected
		rewind=1;Info<<endl;continue;
		} 
	a0=getSelOutput(freak);
	for (j=0; j<nxyz;j++) {for (i=0;i<nsel;i++) { if (freak.spc[i*nxyz+j]<1e10) {species[i*nxyz+j] = freak.spc[i*nxyz+j];} } }
	
	Info << "phreeqc done "<<endl;

	auto finish = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> dt = finish - start;dure = dure+dt.count();std::cout<<"dure"<<dure<<"\n";
	
	// #############  transfer back to C 
	
	forAll(Cw,ic) {
		// if transfered as source term for mobile species, but immobile stay here
		if (immobile[ic]==1) 
		{
			for (j=0;j<nxyz;j++) {Cw[ic]()[ractive[j]] = freak.c[ic*nxyz+j];} // transfer phq -> opf		
		}
		
		//for (j=0;j<nxyz;j++) {if (bcCwi[ractive[j]]==0) {Cw[ic]()[ractive[j]] = freak.c[ic*nxyz+j];} }// transfer phq -> opf		
	}
	
	//Info<<"Cw mid "<<Cw[7]()[int(ncell/2)]<<endl;

	//----------- finding the time step (not for immobile species)
	if (reactionSteps==0) // in this case, we use dtForChem
		{
		dtForChem = 1e12;
		forAll(Cw,ic) // dissolved
			if ((ic>3)&&(immobile[ic]==0))
			{
				dC = 0;dff=0;mxC=0;mnC=1;
				for (j=0;j<nxyz;j++) 
					{mnC=min(mnC,freak.c[ic*nxyz+j]);
					mxC=max(mxC,freak.c[ic*nxyz+j]);}
				for (j=0; j<nxyz;j++)
					if (bcCwi[j]==0) 
					{
					dff = mag(Cw[ic]()[ractive[j]]-freak.c[ic*nxyz+j]);
					dC = max(dC,dff);
					//Info<<"ic, j "<<ic<<" "<<j<<" c "<<Cw[ic]()[ractive[j]]<<" "<<freak.c[ic*nxyz+j]<<" "<<dff<<" "<<dC<<endl;
					}
				dC = dC/(mxC-mnC+SMALL);
				dtForChem = min(dtForChem,dCmax/(max(dC,0)+SMALL)*runTime.deltaTValue()); //Info<<"ic "<<ic<<" dC "<<dC<<" dtForChem "<<dtForChem<<endl;
			} 
		Info<< "dt "<<runTime.deltaTValue()<<" dC "<<dC<<" dtForC "<<dtForC<<" dtForChem " << dtForChem << endl; 
		newDeltaT = min(dtForChem, newDeltaT);
		if ((dtForChem>dtForC*10)&&(tstep>10)) {rSteps=10;} //even if dtForChem is high, keep however a mini value of for security
		else if ((dtForChem>dtForC*2)&&(tstep>10)) {rSteps=static_cast<int>(std::round(dtForChem/dtForC));}
		else {rSteps=1;}
		}
	Info<<"new deltat "<<newDeltaT<<endl;

	//************************ gas, read partial pressures (freak.g in atm) set Cg to fraction (freak.g/Cgtot) and set p to sum of Cg
	
	if (ph_gcomp>0) // multiphase
		{
		float r;
		for (j=0; j<nxyz;j++) //should consider ractive
			{
			//Cgtot = 0; Gmtot = 0;
			//forAll(Cg,i) {freak.g[i*nxyz+j] = max(freak.g[i*nxyz+j],0.); Cgtot += freak.g[i*nxyz+j]; Gmtot += freak.gm[i*nxyz+j];}
			forAll(Cg,i) {Cg[i]()[j] = freak.gm[i*nxyz+j]/gvol[j];} 
			if ((flowType==4)&&(activateEbullition==1)){
				//Vmol[j] = freak.gvol[j]/Gmtot; //1.003 strange but needed to have equilibrated pressure
				//r = freak.p[j]/(p[j]/atmPa);//change of gas pressure
				//sw[j]=1-r*freak.gvol[j]/eps[j];
				p[j] = freak.p[j]*atmPa; //.g are gas pressures, tot is the sum
				}
			//sw[j] = freak.wsat[j];
			//if (j<5) {Info<<"p "<<p[j]<<" sw "<<sw[j]<<" gvol "<<gvol[j]<<" Cg "; for (i=0; i<ph_gcomp;i++){Info<<Cg[i]()[j]<<" ";} Info<<endl;}
			}
		}
	
	// find the variation of wsat from nb moles H2O in gas phase, phreeqc considers a volume of 1 dm3
	if (freak.iGwater>-1) //should consider ractive
		{
		double a1=0.;
		//iw = freak.iGwater; done at start
		for (j=0; j<nxyz;j++)
			{
			a1 = max(freak.gm[iw*nxyz+j],0.) - gm_ph[iw*nxyz+j]; // delta gm water in moles in gvol (equiv of 1L of medium volume)
			//a1 = a1 *   //gm_ph = Cg[i]()*gvol/Vmol
			//if (j<5) {Info<< " gm_ph "<< gm_ph[iw*nxyz+j] << " frk "<< freak.gm[iw*nxyz+j] <<" a1 "<<a1<< endl;}
			sw[j] = max(sw_min[j],sw[j] - a1*.01801/eps[j]); 
			}
		} 
		//nb of moles of H2O(g) transformed in water volume (1 mol 18.01 mL at 25°C)
	//for (j=0;j<3;j++) {Info <<" new sw "<<sw[j]<< endl;}
		
	oldTimeReac = runTime.value()*1;
	reactIndicator = 1; //to know after reaction
} //end activate reaction
