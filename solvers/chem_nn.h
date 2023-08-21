//std::vector<int> cell1; 
std::vector<float> nndata; // TORCH only accepts floats ???? (to vbe validated)
std::vector<float> vnn(nvar,0.);
std::vector<float> vnn1(nvar,0.);
std::vector<double> Cmin {0,0,0,0};
std::vector<double> Cmax {5e-4,1e-3,2e-4,5e-6}; // Dce,Pce, Tce, Vc

if (Cwgnn.trained==1) {
	// reads conc and transfer to vector for nn
	std::ofstream outNNdata_eval(cur_dir/"NNdata_eval.txt");
	nndata.resize(nxyz*nvar);
	for (i=0;i<nxyz;i++) //
	{  
		for (int j=0;j<ph_ncomp-4;j++) { 
			//vnn[j]=static_cast<float>((std::log10(std::max(Cw[j+4]()[i],1e-20))+10)/5);
			vnn[j] = (static_cast<float>(Cw[j+4]()[i] - Cmin[j])/(Cmax[j]-Cmin[j]));
			nndata.push_back(vnn[j]); 
			}
		for (const auto &x : vnn) {outNNdata_eval << x << " ";}
		float dt = static_cast<float>(std::log10(runTime.value()-oldTime)/5);
		nndata.push_back(dt);
		outNNdata_eval<<dt<<"\n";
	}
	Cwgnn.setNd(nxyz);Cwgnn.setData(nndata);Cwgnn.eval();
	// transfer Cwgnn output to C
	int icount=0;
	std::ofstream outNNresu_eval(cur_dir/"NNresu_eval.txt");
	for (i=0;i<nxyz;i++) {
		for (int j=0;j<ph_ncomp-4;j++) {
			double x = static_cast<double>(Cwgnn.output[i][j].item<float>());
			//Cw[j+4]()[i] += 5*std::pow(10,x-10);
			Cw[j+4]()[i] += (x-Cmin[j])*(Cmax[j]-Cmin[j]);
			outNNresu_eval << x <<" ";
			//std::cout<<x<<" ";
			icount +=1;
		}
		icount +=1;
		outNNresu_eval <<"\n";//std::cout<<" "<<std::endl;
	}
}

if (Cwgnn.trained==0) { // case of nn not trained, we store values and use phreeqc to calculate chemistry
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
	
	if ((activateNNchemistry==1)&&(istep>3))
	{
		std::cout<<"start cnn ";
		//put data in the model
		float x;
		// makes an approximate vector to select vairable chemcial composition 
		//dimensionedScalar one = ("one",dimVol/dimMass,1.0);
		//volScalarField Clog = ("Clog",dimless,Cw[4]()*one);Clog=log(Clog); // does not work, strange
		/*std::vector<double>  Clog(nxyz,0.); double mnCl;double mxCl = -35.;
		for (i=0;i<nxyz;i++) { 
			for (int j=4;j<ph_ncomp;j++) { Clog[i] += std::log10(std::max(Cw[j]()[i],1e-7)); }
			//if ((i>150)&&(i<200)) {std::cout<<Clog[i]<<" ";}
			if (Clog[i]<mnCl) {mnCl=Clog[i];}
			if (Clog[i]>mxCl) {mxCl=Clog[i];}
		}
		std::cout<<"mn mx C "<<mnCl<<" "<<mxCl<<std::endl;
		for (i=0;i<nxyz;i++) { 
			Clog[i] = (Clog[i]-mnCl)/(mxCl-mnCl);if ((i>150)&&(i<200)) {std::cout<<Clog[i]<<" ";}
		}*/
		//create a random vecor of indices
		std::vector<int> idx(nxyz);
		std::iota(idx.begin(), idx.end(), 0);		
		std::random_device rd;
		std::mt19937 g(rd());
		std::shuffle(idx.begin(), idx.end(), g);
		
		/*std::vector<int> count(32,0); // will add 16 values from the same category of C/Cmax 
		for (i=0;i<nxyz;i++) 
			{
			i1 = idx[i];int catg = static_cast<int>(Clog[i1]*10);
			if (count[catg]<32) {cell1.push_back(i1);count[catg]+=1;nd +=1;}
			}*/
		//random select where conc changed (rchange)
		
		std::ofstream outNNdata(cur_dir/"NNdata.txt");
		int change;int nd=0;
		std::cout<<" nd "<<nd<<std::endl;
		for (int j=0;j<nxyz;j++) //
		{  
			int j1 = idx[j];change =0;//std::cout<<" vnn "<<j;
			for (int i=0;i<ph_ncomp-4;i++) { 
				//vnn[j]=static_cast<float>((std::log10(std::max(Cw[j+4]()[i1],1e-20))+10)/5);
				//float x = static_cast<float>((Cw[j+4]()[i1] - Cmin[j])/(Cmax[j]-Cmin[j]));
				vnn[i] = static_cast<float>(Cw[i+4]()[j1]); // conc before (not modified)
				vnn1[i] = static_cast<float>(freak.c[(i+4)*nxyz+j1]); // conc after phq
				//std::cout<<" "<<vnn[i]<<" "<<vnn1[i];
				//if ((vnn1[j]-vnn[j])/(vnn[j]+1e-20)>1e-4) {change=1;}
				//nndata.push_back(x); 
			}
			if (rchange[j1]>0) {
				outNNdata << j1 <<" ";
				for (const auto &x : vnn) outNNdata << x << " ";
				for (const auto &x : vnn1) outNNdata << x << " ";
				float dt = static_cast<float>(runTime.value()-oldTime);
				//nndata.push_back(dt);
				//std::cout<<" "<<dt<<std::endl;;
				outNNdata <<dt<<"\n";
			}
			nd +=1 ; if (nd>500) {break;}
		}
		//Cwgnn.setNd(nd);
		//Cwgnn.setData(nndata);
		/*Cwgnn.setTarget(nntarget);
		float rmse = Cwgnn.train(); std::cout << "rmse "<<rmse<<std::endl;
		outNNrmse <<" " << rmse << "\n";
		if (rmse<nn_minr) {Cwgnn.trained = 1;}*/

	} // end of store for activateNN=1
	
	
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
