// * * * * * * * * * * * * * *   steady flow * * * * * * * * * * * * * * //
simpleControl simple(mesh);  // new 11/6/21
double x;
int ilay,nc,nc1;
scalarField resid = h*0;
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
	if (flowType==2) // set the effective thickness (for unconfined flow)
		{
			double minsw=1.;double maxsw =0.;double e=5e-4;
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
					else {sw[nc]=std::max(std::min(x,sw_max),sw_min[nc]);}

					minsw = std::min(minsw,sw[nc]);
					maxsw = std::max(maxsw,sw[nc]);
				}
			}
		Info <<"min max sw "<<minsw<<" "<<maxsw<<endl;
		forAll( mesh.owner() , iface ) 
			{
			oCell = mesh.owner()[iface];odCface = mag(oCf[iface]-oC[oCell]);
			nCell = mesh.neighbour()[iface];ndCface = mag(oCf[iface]-oC[nCell]);
			if (abs(mesh.Sf()[iface].z())<mag(mesh.Sf()[iface])*1e-5)  //vertical faces
				{
				//scalar swm = (sw[oCell]+sw[nCell])/2;
				scalar swm = (sw[oCell]*sw[nCell])/(odCface*sw[oCell]+ndCface*sw[nCell])*(ndCface+odCface);
				Mf[iface] = M0f[iface]*swm;
				}
			/*
			else { //sub-horizontal faces
				if (sw[oCell]<0.8) {Mf[iface]= M0f[iface]*1e-3;}//sw[oCell];} 
				//if (sw[nCell]<0.002) {Mf[iface]= M0f[iface]*sw[iface];}//*sw[nCell];} 
				//scalar swm = min(min(sw[nCell],sw[oCell]),1);
				//Mf[iface]= M0f[iface]*swm;
				}
			*/
			}
		 }
	//Info <<"tolerance "<<Foam::tolerances::solverTolerance<<endl;
	// ***************** solving
	fvScalarMatrix hEqnStd
		(-fvm::laplacian(Mf,h) == b*fvOptions(sw,h) ); //+ fvc::div(phiG) min(mult*nbt,1.)*  b*mult*nbt
	fvOptions.constrain(hEqnStd);
	//hEqnStd.relax();
	//Info << hEqnStd.solver().solve() <<endl; //typeid().name()//hEqnStd.lower() works
	hEqnStd.solve(); //Info << "h 5 " << h[5] << endl;
	//h.relax();
	resid = hEqnStd.residual()()*1;
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


