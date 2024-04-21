// * * * * * * * * * * * * * *   steady flow * * * * * * * * * * * * * * //
simpleControl simple(mesh);  // new 11/6/21
double x;
int ilay,nc,nc1;
//resid.dimensions().reset(dimensionSet(0,1,-1,0,0,0,0));

//Time locTime = Foam::autoPtr<Foam::Time>::New() 	;//runTime.New();
//Foam::Time stdTime(Foam::Time::controlDictName, args); // DONT DO IT it does not converge!!!
while (simple.loop(runTime))
{
Info << "time = " << runTime.timeName() <<  "  deltaT = " <<  runTime.deltaTValue() << endl;

	Uw.correctBoundaryConditions();
	forAll(mesh.boundary(),patchi)
		{
		if (isA< fixedValueFvPatchField<vector> >(Uw.boundaryField()[patchi]))
			{ phiw.boundaryFieldRef()[patchi] = Uw.boundaryField()[patchi] & mesh.Sf().boundaryField()[patchi]; }
		} 
	// ************** set the effective thickness (for unconfined flow)
	if (flowType==2)
		{
			double e=1e-5;
			for (ilay=0;ilay<nlay;ilay++) //ilay will be used from top
			{ 
				for (i=0;i<ncell_lay;i++) 
				{
					nc = (nlay-ilay-1)*ncell_lay + i;
					//if ((ilay>0) && (sw[(nlay-ilay)*ncell_lay + i]>1e-3)) {sw[nc]=1;}
					//else
					//{
					x = (h[nc]-zbot[nc])/thk[nc];
					//if (x<0.02) {sw[nc] = 1e-5+std::exp(x*331-10.5);} 
					//else if (x>0.98) {sw[nc] = 1-1e-5-std::exp((1-x)*331-10.5);}
					//if ((x>1)||(x<0)) {Info<<"x>1 or <0 cell x "<<nc<<" "<<x<<endl;}
					if (x<e && x>0) {sw[nc]=pow(x,2)/(4*e)+x/2+e/4;} // same eq as in mdf USg (transport manual fg 21)
					else if (x>1-e && x<1) {sw[nc]=1-pow((1-x),2)/(4*e)-(1-x)/2-e/4;}
					else {
					sw[nc]=std::max(std::min(x,sw_max),sw_min[nc]);}
				}
			}
		//Info <<"min max sw "<<minsw<<" "<<maxsw<<endl;
		forAll( mesh.owner() , iface ) 
			{
			oCell = mesh.owner()[iface];odCface = mag(oCf[iface]-oC[oCell]);
			nCell = mesh.neighbour()[iface];ndCface = mag(oCf[iface]-oC[nCell]);
			if (abs(mesh.Sf()[iface].z())<mag(mesh.Sf()[iface])*1e-5)  //vertical faces
				{
				//scalar swm = min(sw[oCell],sw[nCell]);
				scalar swm = (sw[oCell]*sw[nCell])/(odCface*sw[oCell]+ndCface*sw[nCell])*(ndCface+odCface);
				Mf[iface] = M0f[iface]*swm;
				}
		else {Mf[iface] = M0f[iface];}
		/*else //other faces if sw=sw_min Mf=Small
			{
			if ((sw[oCell]<sw_min[oCell]*1.1)&&(sw[nCell]<sw_min[nCell]*1.1)) {Mf[iface] = M0f[iface]*sw_min[oCell];}
			}*/
			}
		 }
	//Info <<"tolerance "<<Foam::tolerances::solverTolerance<<endl;
	// -------------- conditions  ---------------------- (can be simplified only the 1st time is needed)
	float time = mesh.time().value(); //Info<<"chw size "<<cellsHw.size()<<" chfx size "<<cellsHfx.size()<<endl;//const Time&
	//get the well discharge
	hSu=hSu*0; hSp=hSp*0;size_t i; //previous cells with values must be destroyed
	//Info << "Hwdata "<<cellsHwData[0]<<" "<<cellsHwData[1]<<" "<<cellsHwData[2]<<" "<<cellsHwData[3]<<endl;
	if (cellsHwel.size()>0)
	{
		std::cout<<" hsuw ";
		for (i=0; i<cellsHwel.size(); i++) 
			{i1=cellsHwel[i]; v=mesh.V()[i1];hSu[i1] = cellsHwelData[i*4+2]/v ;} // discharge must be in s-1
	} Info<<"\nHw read ";
	// set the values in the matrix
	//get the ghb conditions
	if (cellsHghb.size()>0)
	{
		std::cout<<" hsughb ";
		for (i=0; i<cellsHghb.size(); i++) 
			{i1=cellsHghb[i]; v=mesh.V()[i1];hSu[i1] = cellsHghbData[i*4+2]/v ; hSp[i1] = cellsHghbData[i*4+3]/v ; }
	} Info<<"Hghb read ";
	//get the riv conditions
	if (cellsHriv.size()>0)
	{
		for (i=0; i<cellsHriv.size(); i++) 
			{i1=cellsHriv[i];v=mesh.V()[i1];hSu[i1] = cellsHrivData[i*4+2]/v ; hSp[i1] = cellsHrivData[i*4+3]/v ; }
	} Info<<"Hriv read ";
	//get the drn conditions
	if (cellsHdrn.size()>0)
	{
		std::cout<<" hdrn ";
		for (i=0; i<cellsHdrn.size(); i++) 
			{
			i1 = cellsHdrn[i];v=mesh.V()[i1];//Info<<"i1 "<<i1<<" h "<<h()[i1]<<endl;
			hSu[i1] = cellsHdrnData[i*4+2]/v ; hSp[i1] = cellsHdrnData[i*4+3]/v ; 
			//hSp[i1] = cellsHdrnData[2+(icount+i)*4+3] ;
			//if (i<5) {Info<<"drn i1 "<<i1<<" h(i1) "<<h()[i1]<<" hSu "<<hSu[i1]<<" hsp "<<hSp[i1]<<endl;}
			//scalar dff = hSu[i1]-hSp[i1]*h()[i1];
			//if (dff>0) {hSu[i1]/=1000;hSp[i1]/=1000;}
			}
	} Info<<"Hdrn read ";
	// recharge conditions	
	float dv =0;	
	if (cellsHrch.size()>0)
	{
		for (i=0; i<cellsHrch.size(); i++) 
			{
			i1 = cellsHrch[i];
			for (ilay=0;ilay<nlay;ilay++) { if (sw[i1-ilay*ncell_lay]>0.999) {break;} ;} //search for the first confined layer from the top
			ilay1 = max(0, ilay-1); // ilay1 from top
			nc = i1-ilay1*ncell_lay;
			hSu[nc] += cellsHrchData[i*4+2]/mesh.V()[nc] ; // injeciton of recharge in the lowest unconfined layer ) /mag(psi()[cells_[i]])
			} 
	} Info<<"Hdrn read "<<endl;
	//fixed head cells
	/*float cnd0 = 1e-4; float cnd=0;
	if (cellsHfx.size()>0)
	{
		for (int i=0; i<cellsHfx.size(); i++) 
			{i1 = cellsHfx[i]; cnd = cnd0*(mesh.V()[i1]/thk[i1]);
			hSu[i1] = cellsHfxData[2+i*4+2]*cnd; hSp[i1] = cnd; }
			//if (i<5) {Info<<i<<" ifx "<<i1<<" h"<<h()[i]<<" cnd "<<cnd<<" hsu "<<hSu[i1]<<" hsp "<<hSp[i1]<<endl; }}
	} //Info<<"hfx read "<<endl; */


	// ***************** solving
	fvScalarMatrix hEqnStd
		(-fvm::laplacian(Mf,h) - hSu + fvm::Sp(hSp,h) ); //+ fvc::div(phiG) min(mult*nbt,1.)*  b*mult*nbt  == b*fvOptions(sw,h) +fvm::SuSp(hSp,h)
	//fvOptions.constrain(hEqnStd);
	//hEqnStd.relax();
	//get the fixed head conditions
	if (cellsHfix.size()>0)
	{
		for (i=0; i<cellsHfix.size(); i++) {hEqnStd.setReference(cellsHfix[i],cellsHfixData[i*4+2]);}//Info<<"setRef "<<cellsHfx[i]<<endl;}
	}
	//Info << hEqnStd.solver().solve() <<endl; //typeid().name()//hEqnStd.lower() works
	hEqnStd.solve(); //Info << "h 5 " << h[5] << endl;
	//h.relax();
	residu = average(hEqnStd.residual()());
	phiw = hEqnStd.flux() ;
	//phiw = - ((Mf * fvc::interpolate(fvc::grad(h)) ) & mesh.Sf() ); //phiG- 
	Uw = fvc::reconstruct(phiw);
	
			Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
			<< "  ClockTime = " << runTime.elapsedClockTime() << " s"
			<< nl << endl;

}
//std::ofstream hresFile(cur_dir/"hresiduals");
//for (i==0;i<cactive.size();i++) {outFile << CwEqnStd.residual()()[i] << "\n";}
//for (const auto &x : resid) {hresFile << x << "\n";}
//hresFile.close();


