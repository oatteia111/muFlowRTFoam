double x;int ilay,nc,nc1;Info<<"nlay "<<nlay<<" ncl "<<ncell_lay<<endl;//double pi = 3.141592654;
double e=5e-4;
//Info << "iter picard n "<< iterPicard << endl;
h.storePrevIter();Info<<"flowtype "<<flowType<<endl;

if (flowType == 2) // set the effective thickness (for unconfined flow)
	{
	for (ilay=0;ilay<nlay;ilay++) //ilay will be used from top
		{ 
		for (i=0;i<ncell_lay;i++) 
			{
				nc = (nlay-ilay-1)*ncell_lay + i;
				//if ((ilay>0) && (sw[(nlay-ilay)*ncell_lay + i]>1e-3)) {sw[nc]=1;}
				//else
				//{
				x = (h[nc]-zbot[nc])/thk[nc];
				//if (x<0.02) {sw[nc] = 1e-5+std::exp(x*331-10.5);} //1e-4+std::exp(x*150-6.9);} 
				//else if (x>0.98) {sw[nc] = 1-1e-5-std::exp((1-x)*331-10.5);}
				//if (x<0.0005) {sw[nc]=1e-4;}
				//else if (x>0.9995) {sw[nc]=1-1e-4;}
				//else {sw[nc]=x;}
				if (x<e && x>0) {sw[nc]=pow(x,2)/(4*e)+x/2+e/4;} // same eq as in mdf USg (transport manual fg 21)
				else if (x>1-e && x<1) {sw[nc]=1-pow((1-x),2)/(4*e)-(1-x)/2-e/4;}
				else {sw[nc]=std::max(std::min(x,sw_max),sw_min[nc]);}

				//}
				//if (i==2003) {Info<< "cell "<< i <<" lay " << ilay <<" new c "<<nc<< " thk "<<thk[nc]<<" h "<<h[nc]<< " sw "<<sw[nc]<<endl;}

			}
		}
	forAll( mesh.owner() , iface ) 
		{
		// Get the owner and neighbour cell index
		oCell = mesh.owner()[iface];odCface = mag(oCf[iface]-oC[oCell]);
		nCell = mesh.neighbour()[iface];ndCface = mag(oCf[iface]-oC[nCell]);
		if (abs(mesh.Sf()[iface].z())<mag(mesh.Sf()[iface])*1e-5)  //vertical faces
			{
			//scalar swm = (sw[oCell]+sw[nCell])/2;
			//scalar swm = 2*(sw[oCell]*sw[nCell])/(sw[oCell]+sw[nCell]);
			scalar swm = (sw[oCell]*sw[nCell])/(odCface*sw[oCell]+ndCface*sw[nCell])*(ndCface+odCface);
			Mf[iface] = M0f[iface]*swm;
			}
		/*
		else {
			//if (sw[oCell]<0.995) {Mf[iface]= M0f[iface]*1e-3;}//sw[oCell];} 
			if (sw[nCell]<0.002) {Mf[iface]= M0f[iface]*sw[iface];}//*sw[nCell];} 
			//scalar swm = min(min(sw[nCell],sw[oCell]),1);
			//Mf[iface]= M0f[iface]*swm;
			}
		*/
		//Info << " ifac "<<iface << " K "<<Kf[iface] <<endl;
		}
		
	}
forAll(stor(), celli) { 
	if (sw[celli]<0.999) {stor[celli] = eps[celli]/mesh.V()[celli];}// //Info << stor[celli]<<endl;}  /thk[celli]
	else {stor[celli] = stor0/mesh.V()[celli] ;} //   /thk[celli]
	//stor[celli] = (stor0+eps[celli]*(0.49695+1.01/pi*std::atan((0.96-sw[celli])*1000)))/h[celli];// approx formula tomicmic a smoothed step
	};  Info << "stor0 " << stor0 <<" stor(4) " << stor[4] << " eps " << eps[4]<<" thk " <<thk[4]<<" vol "<<mesh.V()[4]<<endl;
stor.dimensions().reset(dimless/dimLength);
// ***************** solving
//Info << "Kf 2026 "<< Kf[2026]<< " " <<Kf[2026]<< endl;

fvScalarMatrix hEqn
	( stor*fvm::ddt(h) - fvm::laplacian(Mf,h) == b*fvOptions(sw,h) ); // + fvc::div(phiG)
fvOptions.constrain(hEqn);
hEqn.solve(); //Info << "h 500 " << h[500] << " iter "<<runTime.solverInfo()<<endl;//nIterations() //hEaqn.solver() ok //SolverPerformance() not ok
Info<<"end solve"<<endl;
phiw =  hEqn.flux() ;
//phiw = phiG- ((Mf * fvc::interpolate(fvc::grad(h)) ) & mesh.Sf() )  ;Info<<"phiw calc"<<endl;
//phiw = phiG - (Mf * fvc::grad(h))* mesh.magSf();Info<<" phiw calc"<<endl;
