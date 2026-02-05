hp.storePrevIter();
sw.storePrevIter();

// -------------- conditions  ----------------------
	hSu=hSu*0; hSp=hSp*0;size_t i; //previous cells with values must be destroyed
	time = mesh.time().value() ; //Info<<"chw size "<<cellsHw.size()<<" chfx size "<<cellsHfx.size()<<endl;//const Time&
	//get the well discharge
	if (cellsHwel.size()>0)
	{
		if (time>tHwel) {
			output = getCbuffer("hwel",itHwel,cellsHwel.size()); 
			tHwel = output.t; cellsHwelData=output.d; itHwel+=1;
			tnext = min(tnext,tHwel);
			}
		for (i=0; i<cellsHwel.size(); i++) {
			i1=cellsHwel[i];qHwel[i] = cellsHwelData[i*4+2] ; hSu[i1] = qHwel[i]/mesh.V()[i1];
			//Info<<"q "<<qHwel[i]<<" "<<mesh.V()[i1]<<endl;
			}
		//Info<<"tHwell "<<tHwel<<" tnext "<<tnext<<endl;
	}
	//----- get the drn conditions
	if (cellsHdrn.size()>0)
	{
		if (time>tHdrn) {
			output = getCbuffer("hdrn",itHdrn,cellsHdrn.size()); 
			tHdrn = output.t; cellsHdrnData=output.d; itHdrn +=1;
			tnext = min(tnext,tHdrn);
			}
		for (i=0; i<cellsHdrn.size(); i++) 
			{
			i1 = cellsHdrn[i]; v=mesh.V()[i1];
			hSu[i1] = cellsHdrnData[i*4+2]/v ; hSp[i1]=cellsHdrnData[i*4+3]/v ;
			if (hp()[i1]< cellsHdrnData[i*4+2]/cellsHdrnData[i*4+3] ) {hSu[i1]/=1e9;hSp[i1]/=1e9;} // almost no flow if head below fixed one
			}
	}  //Info<<"Hdrn read ";
	//---- recharge conditions		
	if (cellsHrch.size()>0)
	{
		if (time>tHrch) {
			output = getCbuffer("hrch",itHrch,cellsHrch.size()); 
			tHrch = output.t; cellsHrchData=output.d; itHrch +=1;
			tnext = min(tnext,tHrch);
			}
		float frech=1;float pi=3.141592654;float hplim = -1.0/lunits; //frech is used to stop recharge or evapo when hp is too high or too low
		for (i=0; i<cellsHrch.size(); i++) 
			{
			i1 = cellsHrch[i];clayRch[i] = i1;qHrch[i] = cellsHrchData[i*4+2] ; // here the injection is at the exact place of recharge (not like in unconfined)
			//Info<<i<<" ic "<<cellsHrch[i]<<" i lay "<<clayRch[i]<<" rch "<<qHrch[i]<<" tHrch "<<tHrch<<" tnext "<<tnext<<endl;
			if (hSu[i1]==0) { // apply recharge only if there is nothing else in the cell
				if (qHrch[i]>0) {frech=1;}
				else 
					{
					if (hp()[i1]< hplim) {frech=0;}//Info<<"frech "<<i1<<" "<<hp()[i1]<<" "<<frech<<endl;} // hp<-3 m stop all evt
					else if (hp()[i1]<0) {frech = 0.5*(1+std::erf((hp()[i1]+0.4/lunits)/20));}//0.5+std::atan((hp()[i1]+0.4/lunits)/8.0)/pi;} //correction to limit head, if cm, lunits= 0.01 -> 40 cm head
					else {frech = 1;}
					}
				hSu[i1] = frech*qHrch[i]/mesh.V()[i1]; 
				//Info<<"rch "<<lunits<<" "<<i1<<" "<<hp()[i1]<<" "<<frech<<endl;
				}
			//Info<<"q rch "<<i<<" " <<qHrch[i]<<endl;
			} 
	}
	//hp fixed
	//std::cout<<"hfix size "<<cellsHfix.size()<<"\n";
	fvScalarMatrix hpEqn
	(stor*fvm::ddt(hp) - fvm::laplacian(Mf,hp) + fvc::div(phiGr) ==  hSu - fvm::SuSp(hSp,hp)); //

	if (cellsHfix.size()>0)
	{
		if (time>=tHfix) {
			output = getCbuffer("hfix",itHfix,cellsHfix.size()); 
			tHfix = output.t; cellsHfixData=output.d; itHfix +=1;
			tnext = min(tnext,tHrch);
			}
		//for (i=0; i<cellsHfix.size(); i++) {i1 = cellsHfix[i];hSu[i1]=-cellsHfixData[i*4+2]*1e3;hSp[i1]=1e3;std::cout<<"h fix "<<i1<<" "<<hSu[i1]<<" "<<hSp[i1];}
		for (i=0; i<cellsHfix.size(); i++) {hpEqn.setReference(cellsHfix[i],cellsHfixData[i*4+2]);}
	}

	//scalar converge = 0;
	
	//if (hp().mesh().equationRelaxationFactor("hp")) {hpEqn.relax(hp().mesh().equationRelaxationFactor("hp"));}
	//hpEqn.solve(); 

	residu = hpEqn.solve().initialResidual(); //initial residual does not seem to include hSu??

//phiw =  phiGr - ((Mf * fvc::interpolate(fvc::grad(hp)) ) & mesh.Sf() ) ;//WRONG
phiw ==  phiGr - (Mf * fvc::snGrad(hp)) * mesh.magSf();//fvc::interpolate(sw) *
Uw = fvc::reconstruct(phiw);
//phiw=hpEqn.flux();
//Info<<"Mf(end) "<<Mf[ncell-1]<<" phiw " <<phiw[ncell-1]<<endl;
//retrieving fluxes at boundaries
volScalarField q1 = fvc::div(phiw);
for (i=0;i<cellsHfix.size();i++) {i1=cellsHfix[i];qHfix[i] = q1[i1]*mesh.V()[i1];} //Info<<" heq hfix "<<i<<" "<<i1<<" "<<qHfix[i]<<endl;
//for (i=0;i<cellsHwel.size();i++) {i1=cellsHwel[i];qHwel[i] = q1[i1]*mesh.V()[i1];} //already done
//for (i=0;i<cellsHrch.size();i++) {i1=cellsHrch[i];qHrch[i] = q1[i1]*mesh.V()[i1];}
for (i=0;i<cellsHdrn.size();i++) {i1=cellsHdrn[i];qHdrn[i]=(hSu[i1]-hp()[i1]*hSp[i1])*mesh.V()[i1];}
//for (i=0;i<cellsHghb.size();i++) {i1=cellsHghb[i];qHghb[i]=(hSu[i1]-hp()[i1]*hSp[i1])*mesh.V()[i1];}

//Info<<"hp -relax";for (int i=0;i<6;i++) {Info<<hp[i*12]<<" ";};Info<<hp[98]<<" "<<hp[99]<<endl;
//hp.relax(); // mix with previous solution to stabilise, relax factor in fvsolu
//if(timeStepControl == "Picard")  {resPicard=residu;}
