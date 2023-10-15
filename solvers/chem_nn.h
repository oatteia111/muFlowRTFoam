double x;
int nmax=750;int k;
std::vector<int> idx(nxyz);
std::vector<float> nndata;
std::vector<float> nntarget;

//##################- case of trained network ###########################--
if ((activateNNchemistry==1)&&(istep>=3)) { //(Cwgnn.trained==1)&&
	//Cwgnn.trained=0;
	// reads conc and transfer to vector for nn
	std::ofstream outNNdata_eval(cur_dir/"NNdata_eval.txt");
	nndata.resize(0); // TORCH only accepts floats ???? (to vbe validated)
	double x;float x1;
	//get selected species (only for 1st tstep)
	std::cout<<"start eval ";
	if (istep==2)
		{
		freak.getSelOutput();	
		S_ph.resize(nxyz*nsel); // S_ph will store all data on surface
		for (j=0;j<nxyz;j++)
			{
			S_ph[j]=freak.spc[j];for (i=1;i<nsel;i++) {S_ph[i*nxyz+j]=freak.spc[(i+1)*nxyz+j];}
			}
		}
	std::cout<<"1st step \n";
	std::cout<<"nxyz "<<nxyz<<" nncols "<<nncols<<" llog0 "<<llog[0]<<" S_ph sz 0 "<<S_ph.size()<<" "<<S_ph[0]<<" ncom "<<ph_ncomp<<"\n";
	//for (i=0;i<ph_ncomp;i++) {std::cout<<Cw[i]().size()<<" ";}
	//rescale data and send to nndata
	for (j=0;j<nxyz;j++)
		{  
		//if (j==0) {std::cout<<"Cw "; for (i=0;i<ph_ncomp;i++) {std::cout<<Cw[i]()[j]<<" ";}}
		//if (j==0) {std::cout<<"Cw ";}
		for (i=0;i<nncols;i++) { 
			if (i<ph_ncomp-4) {x=Cw[i+4]()[j];} //dissolved species
			else {x=S_ph[(i-ph_ncomp+4)*nxyz+j];} //selected out species
			x = std::max(x,1e-12);
			if (llog[i]==1) 
				{x1=static_cast<float>((std::log10(x)-Cmin[i])/(Cmax[i]-Cmin[i]));}
			else 
				{x1=static_cast<float>((x-Cmin[i])/(Cmax[i]-Cmin[i]));}
			if ((j==10)&&(i<ph_ncomp-4)) {std::cout<<" i "<<i<<" cw "<<x<<" x1 "<<x1<<" mx,mn "<<Cmax[i]<<" "<<Cmin[i]<<"\n";}
			outNNdata_eval<<x1<<" ";nndata.push_back(x1);
			}
		outNNdata_eval<<"\n";
		}
	outNNdata_eval.close();
	std::cout<<"in eval \n";
	Cwgnn.setNd(nxyz);Cwgnn.setData(nndata);Cwgnn.eval();
	// transfer Cwgnn output to Cwi and S_ph (we don't need log here, output is dC)
	std::ofstream outNNresu_eval(cur_dir/"NNresu_eval.txt");
	float dt = static_cast<float>(runTime.value()-oldTime); Info<<"dt "<<dt<<endl;
	for (j=0;j<nxyz;j++) {
		for (int i=0;i<nncols;i++) {
			x = static_cast<double>(Cwgnn.output[j][i].item<float>());outNNresu_eval << x <<" "; //trget is dC scaled
			if (i<ph_ncomp-4) {Cw[i+4]()[j] = std::max(Cw[i+4]()[j]+Ymin[i]+ x*(Ymax[i]-Ymin[i]),0.);}
			else {
				S_ph[(i-ph_ncomp+4)*nxyz+j] = std::max(S_ph[(i-ph_ncomp+4)*nxyz+j]+Ymin[i]+ x*(Ymax[i]-Ymin[i]),0.);
				}
			if ((j==10)&&(i<ph_ncomp-4)) {std::cout<<" i "<<i<<" x "<<x<<" cw "<<Cw[i+4]()[j]<<" Y "<<Ymax[i]<<" "<<Ymin[i]<<"\n";}
			//if ((j==10)&&(i>=ph_ncomp-4)) {std::cout<<" i "<<i<<" sph "<<S_ph[(i-ph_ncomp+4)*nxyz+j]<<" Y "<<Ymax[i]<<" "<<Ymin[i]<<"\n";}
		}
		outNNresu_eval <<"\n";//std::cout<<" "<<std::endl;
	}
	outNNresu_eval.close();
}

//########################- NN not trained, PHREEQC ###########################
else
{ // case of nn not trained, we store values and use phreeqc to calculate chemistry
	// find the cells where the chcemistry has changed to calculate there
	//icnt = 0;
	
	//creates rchange where conc has changed since last time step	
	deltaTchem = transportProperties.lookupOrDefault<scalar>("deltaTchem",86400);Info<<"dtchem in reac "<<deltaTchem<<endl;
	for (i=4; i<ph_ncomp; i++)
		{
		for (j=0; j<nxyz;j++)
			{
			if (abs(c_ph[i*nxyz+j]-Cw[i]()[ractive[j]])/(c_ph[i*nxyz+j]+1e-20)>1e-4 && sw[ractive[j]]>sw_min[j]*1.5) {icnt++;rchange[j] = 1.;}
			}
		} 
	for (i=0; i<ph_gcomp; i++)
		{
		for (j=0; j<nxyz;j++)
			{
			if (abs(gm_ph[i*nxyz+j]-Cg[i]()[ractive[j]])/(gm_ph[i*nxyz+j]+1e-20)>1e-5 && sw[ractive[j]]>sw_min[j]) {rchange[j] = max(rchange[j],1.);}
			}
		} 
	// transfer to phreeqc for c and g, for g we send moles and not pressures, phreeqc does not know the cell volume it considers to be one
	//int icnt = 0;
	for (j=0; j<nxyz;j++) { forAll(Cw,i) {c_ph[i*nxyz+j] = Cw[i]()[ractive[j]];} } 
	if (flowType == 4)
		{
		for (j=0; j<nxyz;j++) { 
			forAll(Cg,i) {
			gm_ph[i*nxyz+j] = max(Cg[i]()[j]*gvol[j]/phreeqcVm[j],1e-6);//Cg in fraction and Vm in L/mol
				}	
			} 
		}
	else 
		{
		for (j=0; j<nxyz;j++) { 
			forAll(Cg,i) {
				gm_ph[i*nxyz+j] = max(Cg[i]()[j]*gvol[j]/phreeqcVm[j],1e-16);
				} //Cg in fraction and Vm in L/mol
			} 
		}	
	Info<<" phqVm[0] "<<phreeqcVm[0]<<" "<<endl;
	
	//auto start = std::chrono::high_resolution_clock::now();
	//################# RUN PHREEQC   ################
	// set saturations using rchange
	for (j=0; j<nxyz;j++) {rchange[j] *= sw[j];}//Info<<"rch "<<rchange[j]<<endl;}
	freak.setGvol(gvol); // set gas volume in phreeqc
	freak.setWsat(rchange); // rchange for the calculation doamin, with 0 outside, sw saturation
	freak.setC(c_ph);//transfer c_ph to freak : it does not work to send directly to freak.c
	freak.setGm(gm_ph);//transfer gm_ph to freak
	//freak.setP(p_ph);//transfer pressure to freak
	freak.setTstep(runTime.value()-oldTime); //Info<<" this tme "<< runTime.value()<<" old "<<oldTime<<endl;//the calculation time shall include all time since las phreeqc run
	Info << "running phreeqc dt "<<runTime.value()-oldTime<<endl;
	freak.run();
	freak.getSelOutput();
	Info << "phreeqc done "<<endl;
		
	// write to intermediate file, input (for gases)
	if (ph_gcomp>0) {
		std::ofstream inPhq(cur_dir/"phq_input.txt");
		for (j=0; j<nxyz;j++) { inPhq<<j<<" "<<p[j]<<" "<<rchange[j]<<" "<<gvol[j]<<" "; for (i=0; i<ph_gcomp;i++) {inPhq<< gm_ph[i*nxyz+j] <<" ";} inPhq<<"\n"; }
		}
	
	//auto finish = std::chrono::high_resolution_clock::now();
	//std::chrono::duration<double> dt = finish - start;dure = dure+dt.count();
	
	//###############- get conc and send them to data and target  #######################################-
	if ((activateNNchemistry==1)&&(istep<3))
	{	
		std::cout<<"start cnn ";
		//create a random vecor of indices
		//std::iota(idx.begin(), idx.end(), 0);		
		//std::shuffle(idx.begin(), idx.end(), shuf);
		
		std::ofstream outNNdata(cur_dir/"NNdata.txt");
		//recalculate Cmin and Cmax (considering lo)
		for (j=0;j<nxyz;j++)
			{
			for (i=0;i<nncols;i++)
				{
				if (i<ph_ncomp-4) {x = c_ph[i*nxyz+j];}
				else {x=S_ph[i*nxyz+j];}
				x = std::max(x,1e-12);
				if (llog[i]==1) {Cmin[i]=std::min(Cmin[i],std::log10(x));Cmax[i]=std::max(Cmax[i],std::log10(x));}
				else {Cmin[i]=std::min(Cmin[i],x);Cmax[i]=std::max(Cmax[i],x);}
				}
			}
		std::cout<<"Cmin "; for (i=0;i<nncols;i++) {std::cout<<Cmin[i]<<" ";} ; std::cout<<"\n";
		std::cout<<"Cmax "; for (i=0;i<nncols;i++) {std::cout<<Cmax[i]<<" ";} ; std::cout<<"\n";
		//transform the variables
		std::vector<float> data(nxyz*nncols);
		for (j=0;j<nxyz;j++)
			{
			for (i=0;i<nncols;i++) 
				{
				if (i<ph_ncomp-4) {x=c_ph[(i+4)*nxyz+j];}
				else {x=S_ph[(i-ph_ncomp+4)*nxyz+j];}
				x = std::max(x,1e-12);
				if (llog[i]==1) {x=static_cast<float>((std::log10(x+cdff[i]/1e12)-Cmin[i])/(Cmax[i]-Cmin[i]));}
				else {x=static_cast<float>((x-Cmin[i])/(Cmax[i]-Cmin[i]));}
				data[i*nxyz+j]=x;
				}
			}
		std::ofstream outNNtarget(cur_dir/"NNtarget.txt");
		freak.getSelOutput();std::vector<double>s_dff(nxyz*nncols); //base
		for (j=0;j<nxyz;j++)
			{
			for (i=0;i<ph_ncomp-4;i++) {s_dff[i*nxyz+j] = freak.c[(i+4)*nxyz+j]-c_ph[(i+4)*nxyz+j];}
			i1=ph_ncomp-4;s_dff[i1*nxyz+j]=freak.spc[j]-S_ph[j];
			for (i=1;i<nsel;i++) {i1=ph_ncomp-4+i;s_dff[i1*nxyz+j]=freak.spc[(i+1)*nxyz+j]-S_ph[i*nxyz+j];}
			}
		for (j=0;j<nxyz;j++)
			{
			for (i=0;i<nncols;i++) {x = s_dff[i*nxyz+j];Ymin[i]=std::min(Ymin[i],x);Ymax[i]=std::max(Ymax[i],x);}
			}
		//we need to find the indices where the concentrations change linked to reaction are important
		//test reorder
		std::vector<int> id0={0,1,2,3,4,5};
		std::vector<double> y0={44.1,12.,3.,28.5,52.1,0.78};
		sort( id0.begin(),id0.end(), [&](int i,int j){return y0[i]>y0[j];} );
		std::cout<<"test ";for (i=0;i<6;i++) {std::cout<<id0[i]<<" ";}std::cout<<"\n";
		std::vector<double> Ysum(nxyz);
		for (j=0;j<nxyz;j++) {for (i=0;i<nncols;i++) {Ysum[j] += s_dff[i*nxyz+j];} }
		std::iota(idx.begin(), idx.end(), 0);	
		sort( idx.begin(),idx.end(), [&](int i,int j){return Ysum[i]>Ysum[j];} ); //sorted from higher to lower_bound

		//then write these data in nndata
		for (k=0;k<nmax;k++)
			{
			j = idx[k];
			for (i=0;i<nncols;i++) {x=data[i*nxyz+j];outNNdata<<x<<" ";nndata.push_back(x);}
			outNNdata<<"\n";
			}
			
		for (i=0;i<nncols;i++) {cdff[i]=Ymax[i]-Ymin[i];}
		std::cout<<"Ymin "; for (i=0;i<nncols;i++) {std::cout<<Ymin[i]<<" ";} ; std::cout<<"\n";
		std::cout<<"Ymax "; for (i=0;i<nncols;i++) {std::cout<<Ymax[i]<<" ";} ; std::cout<<"\n";
		//now write the conc differences in target (for 1st phase)
		for (k=0;k<nmax;k++)
			{
			j=  idx[k];
			//for (i=0;i<nncols;i++) {outNNtmpdC<<s_dff[i*nxyz+j]<<" ";}
			for (i=0;i<nncols;i++) {x=static_cast<float>((s_dff[i*nxyz+j]-Ymin[i])/(Ymax[i]-Ymin[i]+cdff[i]/1e12));outNNtarget<<x<<" ";nntarget.push_back(x);}
			outNNtarget<<"\n";//outNNtmpdC<<"\n";
			}

		Cwgnn.setNd(nmax);//Cwgnn.setRunParms({nn_epoc,nn_batch,nn_lr,0.75});
		Cwgnn.setData(nndata);
		Cwgnn.setTarget(nntarget);
		float rmse = Cwgnn.train();std::cout << "rmse "<<rmse<<std::endl;
	} // end of store for activateNN=1
	
	
	//################## end of NN , store phreeqc results ##############################
	// transfer back to C but before keep the previous values for outside domain of calculation
	forAll(Cw,i) {Cw[i]() = Cw[i]().prevIter();}
	if (ph_gcomp>0) {forAll(Cg,i) {Cg[i]() = Cg[i]().prevIter();} }
	forAll(Cw,i) // dissolved
		{
			for (j=0; j<nxyz;j++)
				{
				Cw[i]()[ractive[j]] = freak.c[i*nxyz+j];
				//if (j==imin) {Info<<"ic "<<i<<" imin "<<imin<<" c "<<Cw[i]()[imin]<<endl;}
				}
		} 
	
	// gas, read partial pressures (freak.g in atm) set Cg to fraction (freak.g/Cgtot) and set p to sum of Cg
	if (flowType == 4) // multiphase
		{
		for (j=0; j<nxyz;j++) //should consider ractive
			{
			Cgtot = 0; Gmtot = 0;
			forAll(Cg,i) {freak.g[i*nxyz+j] = max(freak.g[i*nxyz+j],0.); Cgtot += freak.g[i*nxyz+j]; Gmtot += freak.gm[i*nxyz+j];}
			forAll(Cg,i) {Cg[i]()[j] = freak.g[i*nxyz+j]/Cgtot;} // /Cgtot/phreeqcVm;}
			phreeqcVm[j] = gvol[j]/Gmtot;
			p[j] = Cgtot*atmPa;
			//sw[j] = freak.wsat[j];
			}
		for (int i=0; i<ph_gcomp;i++){for (j=0;j<3;j++) {Info <<"g_spc "<< i <<" Cg "<< Cg[i]()[j] <<" sw "<<sw[j]<< endl;}}
		}
	else if (ph_gcomp>0) //unsaturated (not diffrent from above now)
		{
		for (j=0; j<nxyz;j++) //should consider ractive
			{
			Cgtot = 0;
			//phreeqcVm[j] = freak.g[j]*gvol[j]/freak.gm[j];
			forAll(Cg,i) {freak.g[i*nxyz+j] = max(freak.g[i*nxyz+j],0.); Cgtot += freak.g[i*nxyz+j];}
			forAll(Cg,i) {Cg[i]()[j] = freak.g[i*nxyz+j]/Cgtot;} // /Cgtot or Cgtot/phreeqcVm;}
			}
		}
	//Info<<"gvol 20 "<<gvol[20]<<" cg 0 19 "<<Cg[0]()[19]<<" cg 0 20 "<<Cg[0]()[20]<<" cg 0 21 "<<Cg[0]()[21]<<endl;
	// write to phq output file			
	if (ph_gcomp>0) {
		std::ofstream outPhq(cur_dir/"phq_output.txt");
		for (j=0; j<nxyz;j++) { outPhq<<j<<" "<<p[j]<<" "<<freak.gvol[j]<<" "<<phreeqcVm[j]<<" "<<freak.p[j]<<" "<<freak.wsat[j]<<" "; 
			for (i=0; i<ph_gcomp;i++) {outPhq<< freak.g[i*nxyz+j] <<" ";} 
			for (i=0; i<ph_gcomp;i++) {outPhq<< freak.gm[i*nxyz+j] << " ";} 
			outPhq<<"\n"; 
			}
		}

	// find the variation of wsat from nb moles H2O in gas phase, phreeqc considers a volume of 1 dm3
	if (freak.iGwater>-1) //should consider ractive
		{
		double a1=0.;
		//iw = freak.iGwater; done at start
		for (j=0; j<nxyz;j++)
			{
			a1 = max(freak.gm[iw*nxyz+j],0.) - gm_ph[iw*nxyz+j]; // delta gm water 
			//a1 = a1 *   //gm_ph = Cg[i]()*gvol/phreeqcVm
			if (j<5) {Info<< " gm_ph "<< gm_ph[iw*nxyz+j] << " frk "<< freak.gm[iw*nxyz+j] <<" a1 "<<a1<< endl;}
			sw[j] = max(sw_min[j],sw[j] - a1*.01801/eps[j]); 
			}
		} 
		//nb of moles of H2O(g) transformed in water volume (1 mol 18.01 mL at 25°C)
	for (j=0;j<3;j++) {Info <<" new sw "<<sw[j]<< endl;}
}
