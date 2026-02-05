// * * * * * * * * * * * * * *   steady flow * * * * * * * * * * * * * * //
simpleControl simple(mesh);  // new 11/6/21
double x;double e=1e-5;
int ilay,nc,nc1;
//resid.dimensions().reset(dimensionSet(0,1,-1,0,0,0,0));

//Time locTime = Foam::autoPtr<Foam::Time>::New() 	;//runTime.New();
//Foam::Time stdTime(Foam::Time::controlDictName, args); // DONT DO IT it does not converge!!!
residu = 1e5;
while (simple.loop(runTime))  //(abs(residu)>1e-6) 
{
Info << "time = " << runTime.timeName() <<  "  deltaT = " <<  runTime.deltaTValue() << endl;

	/* if they are bc in Uw :-> removed from now one
	Uw.correctBoundaryConditions();
	forAll(mesh.boundary(),patchi)
		{
		if (isA< fixedValueFvPatchField<vector> >(Uw.boundaryField()[patchi]))
			{ phiw.boundaryFieldRef()[patchi] = Uw.boundaryField()[patchi] & mesh.Sf().boundaryField()[patchi]; }
		}
	*/
	// ************** set the effective thickness (for unconfined flow)
	if (flowType==2)
		{
			for (ilay=0;ilay<nlay;ilay++) //ilay will be used from top
			{ 
				for (i=0;i<ncell_lay;i++) 
				{
					nc = ilay*ncell_lay + i;
					x = (h[nc]-zbot[nc])/thk[nc];
					if (x<0) {dryCells.append(nc);}
					if (x<e && x>0) {sw[nc]=pow(x,2)/(4*e)+x/2+e/4;} // same eq as in mdf USg (transport manual fg 21)
					else if (x>1-e && x<1) {sw[nc]=1-pow((1-x),2)/(4*e)-(1-x)/2-e/4;}
					else {sw[nc]=std::max(std::min(x,sw_max[nc]),sw_min[nc]);}
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
		for (i=0; i<cellsHwel.size(); i++) 
			{i1=cellsHwel[i]; v=mesh.V()[i1];hSu[i1] = cellsHwelData[i*4+2]/v ;}
	} //Info<<"\nHw read ";
	// set the values in the matrix
	//get the ghb conditions
	if (cellsHghb.size()>0)
	{
		for (i=0; i<cellsHghb.size(); i++) 
			{i1=cellsHghb[i]; v=mesh.V()[i1];hSu[i1] = cellsHghbData[i*4+2]/v ; hSp[i1] = cellsHghbData[i*4+3]/v ; }
	} //Info<<"Hghb read ";
	//get the riv conditions
	if (cellsHriv.size()>0)
	{
		for (i=0; i<cellsHriv.size(); i++) 
			{i1=cellsHriv[i];v=mesh.V()[i1];hSu[i1] = cellsHrivData[i*4+2]/v ; hSp[i1] = cellsHrivData[i*4+3]/v ; }
	} //Info<<"Hriv read ";
	//get the drn conditions
	if (cellsHdrn.size()>0)
	{
		for (i=0; i<cellsHdrn.size(); i++) 
			{
			i1 = cellsHdrn[i];v=mesh.V()[i1];//Info<<"i1 "<<i1<<" h "<<h()[i1]<<endl;
			hSu[i1] = cellsHdrnData[i*4+2]/v ; hSp[i1] = cellsHdrnData[i*4+3]/v ; 
			if (h()[i1]< cellsHdrnData[i*4+2]/cellsHdrnData[i*4+3] ) {hSu[i1]/=1000.;hSp[i1]/=1000.;} // no flow if head below fixed one
			}
	} //Info<<"Hdrn read ";
	// recharge conditions	
	float dv =0;	
	if (cellsHrch.size()>0)
	{
		for (i=0; i<cellsHrch.size(); i++) 
			{
			i1 = cellsHrch[i];
			for (ilay=0;ilay<nlay;ilay++) {nc=i1-ilay*ncell_lay; if (sw[nc]>sw_max[nc]) {break;} ;} //search for the first confined layer from the top
			ilay1 = max(0, ilay-1); // ilay1 from top
			nc = i1-ilay1*ncell_lay;
			if (hSu[nc]==0) {
				hSu[nc] += cellsHrchData[i*4+2]/mesh.V()[nc] ; // injeciton of recharge in the lowest unconfined layer ) /mag(psi()[cells_[i]])
				}
			} 
	} //Info<<"Hrch read "<<endl;

	// ***************** solving
	fvScalarMatrix hEqnStd
		(-fvm::laplacian(Mf,h) + fvm::Sp(hSp,h) == hSu); //+ fvc::div(phiG) min(mult*nbt,1.)*  b*mult*nbt  == b*fvOptions(sw,h) +fvm::SuSp(hSp,h)
	//fvOptions.constrain(hEqnStd);
	//get the fixed head conditions
	if (cellsHfix.size()>0)
	{
		for (i=0; i<cellsHfix.size(); i++) {hEqnStd.setReference(cellsHfix[i],cellsHfixData[i*4+2]);}//Info<<"setRef "<<cellsHfx[i]<<endl;}
	}

	residu = hEqnStd.solve().initialResidual(); //Info << "h 5 " << h[5] << endl;
	//h.relax();
	//residu = residu + average(hEqnStd.residual()());
	Info<<"myresidu "<<residu<<endl;
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


