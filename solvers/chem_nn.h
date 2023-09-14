//std::vector<int> cell1; 
std::vector<float> vnn(nvar,0.);
std::vector<float> vnn1(nvar,0.);
	std::vector<float> nndata; // TORCH only accepts floats ???? (to vbe validated)
	std::vector<float> nntarget;

double crel;

if ((activateNNchemistry==1)&&(istep>=20)) { //(Cwgnn.trained==1)&&
	//Cwgnn.trained=0;
	// reads conc and transfer to vector for nn
	std::ofstream outNNdata_eval(cur_dir/"NNdata_eval.txt");
	std::cout<<"nxyz "<<nxyz<<"\n";
	std::vector<float> nndata; // TORCH only accepts floats ???? (to vbe validated)
	std::vector<float> nntarget;float x;
	for (j=0;j<nxyz;j++) //
	{  
		for (int i=0;i<ph_ncomp-4;i++) { 
			//vnn[j]=static_cast<float>((std::log10(std::max(Cw[j+4]()[i],1e-20))+10)/5);
			x= static_cast<float>((Cw[i+4]()[j] - Cmin[i])/(Cmax[i]-Cmin[i]));
			nndata.push_back(x);outNNdata_eval << x << " "; //data is Crelative
			}
		outNNdata_eval<<"\n";
	}
	std::cout<<"in eval \n";
	Cwgnn.setNd(nxyz);Cwgnn.setData(nndata);Cwgnn.eval();
	// transfer Cwgnn output to C
	std::ofstream outNNresu_eval(cur_dir/"NNresu_eval.txt");
	float dt = static_cast<float>(runTime.value()-oldTime); Info<<"dt "<<dt<<endl; 
	for (j=0;j<nxyz;j++) {
		for (int i=0;i<ph_ncomp-4;i++) {
			double x = static_cast<double>(Cwgnn.output[j][i].item<float>()); // the target was x=dC/dt*1e9
			Cw[i+4]()[j] = max(Cw[i+4]()[j]+ x*dt/5e9,0.); //model produces a dC (not dCrelative)
			outNNresu_eval << x <<" ";
		}
		outNNresu_eval <<"\n";//std::cout<<" "<<std::endl;
	}
}

if ((activateNNchemistry==0)||(istep<20)) { // case of nn not trained, we store values and use phreeqc to calculate chemistry
	// find the cells where the chcemistry has changed to calculate there
	//icnt = 0;
	// first stores C values to train NN
	
	// -------------now start phreeqc ------------------------
	deltaTchem = transportProperties.lookupOrDefault<scalar>("deltaTchem",86400);Info<<"dtchem in reac "<<deltaTchem<<endl;
	//creates rchange where conc has changed since last time step
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
	//---------------- send data to NN  ----------------------------------------
	/*
	if ((activateNNchemistry==1)&&(istep>3))
	{
		std::cout<<"start cnn ";
		//put data in the model
		float x;
		// makes an approximate vector to select vairable chemcial composition 
		//create a random vecor of indices
		std::vector<int> idx(nxyz);
		std::iota(idx.begin(), idx.end(), 0);		
		//std::random_device rd;
		//std::mt19937 g(rd());
		std::shuffle(idx.begin(), idx.end(), shuf);
		
		if (nstack==0) { //when nstack=0 we reset all files and vectors to 0
			std::ofstream outNNdata(cur_dir/"NNdata.txt");
			std::ofstream outNNtarget(cur_dir/"NNtarget.txt");
			std::vector<float> nndata; // TORCH only accepts floats ???? (to vbe validated)
			std::vector<float> nntarget;
			} //new files

		float sC;int nC0=0;int nd0=0;  // I want 80% of conc at places where conc !=0
		float dt = static_cast<float>(runTime.value()-oldTime);
		for (int j=0;j<nxyz;j++) //
		{  
			int j1 = idx[j];sC =0;//std::cout<<" vnn "<<j;
			for (int i=0;i<ph_ncomp-4;i++) { 
				vnn[i] = static_cast<float>(Cw[i+4]()[j1]); // conc before (not modified)
				sC += vnn[i];
				vnn1[i] = static_cast<float>(freak.c[(i+4)*nxyz+j1]-vnn[i]); // dC diff conc
				if (sC<1e-7) {nC0+=1;}
			} 
			if ((nC0<50)||(sC>1e-7)) { //50 values close to 0
				for (int i=0;i<ph_ncomp-4;i++) { 
					x = (vnn[i] - Cmin[i])/(Cmax[i]-Cmin[i]); // data is relative conc
					nndata.push_back(x);outNNdata << x << " ";
					x = vnn1[i]/dt*5e9; // target is dC/dt*5e9
					nntarget.push_back(x);outNNtarget << x << " ";
				}
				outNNdata<<"\n";outNNtarget<<"\n";
				nd +=1;nd0+=1;
			}
			if (nd0>500) {break;}
		}
		std::cout<<" nd "<<nd<<" nstack "<<nstack<<std::endl;
		if (nstack>9) {
			Cwgnn.setNd(nd);Cwgnn.setRunParms({nn_epoc,nn_batch,nn_lr,0.75});
			Cwgnn.setData(nndata);
			Cwgnn.setTarget(nntarget);
			float rmse = Cwgnn.train(); std::cout << "rmse "<<rmse<<std::endl;
			outNNrmse <<" " << rmse << "\n";
			if (rmse<nn_minr) {Cwgnn.trained = 1;}
			nstack = 0;nd=0;
		}
	nstack += 1;
	} // end of store for activateNN=1
	*/
	
	//------------------ end of NN , store phreeqc results ------------------------------
	// transfer back to C but before keep the previous values for outside domain of calculation
	forAll(Cw,i) {Cw[i]() = Cw[i]().prevIter();}
	if (ph_gcomp>0) {forAll(Cg,i) {Cg[i]() = Cg[i]().prevIter();} }
	forAll(Cw,i) // dissolved
		{
			for (j=0; j<nxyz;j++)
				{
				Cw[i]()[ractive[j]] = freak.c[i*nxyz+j];
				if (j==imin) {Info<<"ic "<<i<<" imin "<<imin<<" c "<<Cw[i]()[imin]<<endl;}
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
