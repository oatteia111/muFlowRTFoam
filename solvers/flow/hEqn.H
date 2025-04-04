/*
Here we calculate the head equation for different conditions (confined, unconfined) the unsaturated is calculated in hpEqn,
using the boundary conditions
the fluxes at boundaries are calculated in the budget file
*/
double x;int ilay,nc,nc1;//Info<<"nlay "<<nlay<<" ncl "<<ncell_lay<<endl;//double pi = 3.141592654;
double e=1e-5;
//Info << "iter picard n "<< iterPicard << endl;
h.storePrevIter();//Info<<"flowtype "<<flowType<<endl;

//for (i=0;i<150;i++) {Info << "cID " <<i<<" "<< mesh.C()[i] << endl;}

// get the multiplying factors from plugin_H 
volScalarField factM = plugH.corrM(mesh,freak);//Info<<"factM(0-2) "<<factM[0]<<" "<<factM[1] <<" "<<factM[2]<<endl;//" Mf "<<Mf<<" M0f "<<M0f<<endl;
surfaceScalarField factMf = fvc::interpolate(factM,"factM");

// -----------------  set the effective thickness for unconfined flow -------------------
if (flowType == 2) 
	{
	for (ilay=0;ilay<nlay;ilay++) //find sw
		{ 
		for (i=0;i<ncell_lay;i++) 
			{
				nc = (nlay-ilay-1)*ncell_lay + i; // top layer 1st
				x = (h[nc]-zbot[nc])/thk[nc];
				if (x<e && x>0) {sw[nc]=pow(x,2)/(4*e)+x/2+e/4;} // same eq as in mdf USg (transport manual fg 21)
				else if (x>1-e && x<1) {sw[nc]=1-pow((1-x),2)/(4*e)-(1-x)/2-e/4;}
				else {
				sw[nc]=std::max(std::min(x,sw_max[nc]),sw_min[nc]);
				}
			}
		}
		
	forAll( mesh.owner() , iface ) //calculate Mf
		{
		oCell = mesh.owner()[iface];odCface = mag(oCf[iface]-oC[oCell]);
		nCell = mesh.neighbour()[iface];ndCface = mag(oCf[iface]-oC[nCell]);
		if (abs(mesh.Sf()[iface].z())<mag(mesh.Sf()[iface])*1e-5)  //vertical faces
			{
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
else 
	{Mf = M0f*1;} //this is done because if not Mf is x each time by factM below

//Info << "stor0 " << stor0 <<" stor(4) " << stor <<" stU "<<storU<<" stC "<<storC<< " eps " << eps[4]<<" thk " <<thk[4]<<" vol "<<mesh.V()[4]<<endl;
forAll(stor(), celli) 
	{ 
	if (sw[celli]<sw_max[nc]) {stor[celli] = storU[celli]/thk[celli];}// //Info << stor[celli]<<endl;}  
	else {stor[celli] = storC[celli]/thk[celli];} //   /thk[celli]
	//stor[celli] = (stor0+eps[celli]*(0.49695+1.01/pi*std::atan((0.96-sw[celli])*1000)))/h[celli];// approx formula tomicmic a smoothed step
	};  //Info << "stor0 " << stor0 <<" stor(4) " << stor[4] << " eps " << eps[4]<<" thk " <<thk[4]<<" vol "<<mesh.V()[4]<<endl;
stor.dimensions().reset(dimless/dimLength);

// -################# - conditions  ---####################---
	hSu=hSu*0; hSp=hSp*0; size_t i;//previous cells with values must be destroyed
	time = mesh.time().value() ;float f; //Info<<"chw size "<<cellsHw.size()<<" chfx size "<<cellsHfx.size()<<endl;//const Time&
	//get the well discharge/recharge
	if (cellsHwel.size()>0)
	{
		if (time>tHwel) {
			output = getCbuffer("hwel",itHwel,cellsHwel.size()); 
			tHwel = output.t; cellsHwelData=output.d; itHwel+=1;
			tnext = min(tnext,tHwel);
			}
		for (i=0; i<cellsHwel.size(); i++)
			{v=mesh.V()[cellsHwel[i]];qHwel[i] = cellsHwelData[i*4+2] ; hSu[cellsHwel[i]] = qHwel[i]/v;}
			//Info<<i<<" ic well "<<cellsHwel[i]<<" qwel "<<qHwel[i]<<endl; }
	} //Info<<"Hfw read ";
	//----- get the ghb conditions
	//int ic0=10;
	if (cellsHghb.size()>0)
	{
		if (time>tHghb) {
			output = getCbuffer("hghb",itHghb,cellsHghb.size()); 
			tHghb = output.t; cellsHghbData=output.d; itHghb +=1;
			tnext = min(tnext,tHghb);
			}
		for (i=0; i<cellsHghb.size(); i++) 
			{i1=cellsHghb[i];v=mesh.V()[i1];hSu[i1] = cellsHghbData[i*4+2]/v ; hSp[i1] = cellsHghbData[i*4+3]/v ;}
	} //Info<<"Hghb read hSu ";//<<hSu[ncell_lay*3+ic0]<<" ";
	//----- get the riv conditions
	if (cellsHriv.size()>0)
	{
		if (time>tHriv) {
			output = getCbuffer("hriv",itHriv,cellsHriv.size()); 
			tHriv = output.t; cellsHrivData=output.d; itHriv +=1;
			tnext = min(tnext,tHriv);
			}
		for (i=0; i<cellsHriv.size(); i++) 
			{i1=cellsHriv[i];v=mesh.V()[i1];hSu[i1] = cellsHrivData[i*4+2]/v ; hSp[i1] = cellsHrivData[i*4+3]/v ;}
			//Info<<" i "<<i<<" i1 "<<i1<<" hsu "<<hSu[i1]<<" hsp "<<hSp[i1]<<endl;}
	} //Info<<"Hriv read "<<endl;
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
			i1 = cellsHdrn[i]; v=mesh.V()[i1]; // missing condition!!!
			hSu[i1] = cellsHdrnData[i*4+2]/v ; hSp[i1]=cellsHdrnData[i*4+3]/v ;
			if (h()[i1]< cellsHdrnData[i*4+2]/cellsHdrnData[i*4+3] ) {hSu[i1]/=1e6;hSp[i1]/=1e6;} // no flow if head below fixed one
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
		for (i=0; i<cellsHrch.size(); i++) 
			{
			i1 = cellsHrch[i];//Info<<" i "<<i<<" i1 "<<i1<<" nla "<<nlay<<" ncl "<<ncell_lay;
			for (ilay=0;ilay<nlay;ilay++) { nc=i1-ilay*ncell_lay; if (sw[nc]>=sw_max[nc]) {break;} ;} //search for the first confined layer from the top
			ilay1 = max(0, ilay-1);// ilay1 from top
			nc = i1-ilay1*ncell_lay;clayRch[i] = nc;qHrch[i] = cellsHrchData[i*4+2] ;
			//if(hSu[nc]!=0) {Info<<"rch i i1 "<<i<<" "<<i1<<" nc "<<nc<<" h(nc) "<<h()[nc]<<" sw[nc] "<<sw[nc]<<" sw_min "<<sw_min[nc]<<" data "<<cellsHrchData[2+(icount+i)*4+2]<<endl;}
			if (hSu[nc]==0) {
				hSu[nc] = qHrch[i]/mesh.V()[nc]; 
				} // injeciton of recharge in the lowest unconfined layer only if no fix cond here
			} 	
	} //Info<<"Hrch read ";//<<" Hghb "<<hSu[ncell_lay*3+ic0]<<" "<<endl;

// ***************** solving
// linked to coupling
Mf *=factMf; 
Info<<"Mf(0-2) "<<Mf[0]<<" "<<Mf[1] <<" "<<Mf[2]<<endl;
fvScalarMatrix hEqn
	( stor*fvm::ddt(h) - fvm::laplacian(Mf,h) == hSu - fvm::SuSp(hSp,h) ); // + fvc::div(phiG) == b*fvOptions(sw,h) 
//fixed head cells
if (cellsHfix.size()>0)
{
	if (time>tHfix) {
		output = getCbuffer("hfix",itHfix,cellsHfix.size()); 
		tHfix = output.t; cellsHfixData=output.d; itHfix +=1;
		tnext = min(tnext,tHfix);
		}
	//std::cout<<"hfix ";
	for (i=0; i<cellsHfix.size(); i++) {hEqn.setReference(cellsHfix[i],cellsHfixData[i*4+2]);}//std::cout<<cellsHfix[i]<<" "<<cellsHfixData[i*4+2]<<" ";}
}
//resid = average(hEqn.residual());

residu = hEqn.solve().initialResidual(); 
Info<<"h(0-2) "<<h[0]<<" "<<h[1] <<" "<<h[2]<<endl;
//phiw= -((Mf*gradHf) & mesh.Sf());
phiw == (-Mf * fvc::snGrad(h)) * mesh.magSf(); //similar to the above one
Uw = fvc::reconstruct(phiw);
volScalarField q1 = fvc::div(phiw);
for (i=0;i<cellsHfix.size();i++) {i1=cellsHfix[i];qHfix[i] = q1[i1]*mesh.V()[i1];}//Info<<" heq hfix "<<i<<" "<<i1<<" "<<qHfix[i]<<endl;}
//for (i=0;i<cellsHwel.size();i++) {i1=cellsHwel[i];qHwel[i] = q1[i1]*mesh.V()[i1];} //already calculated
//for (i=0;i<cellsHrch.size();i++) {i1=cellsHrch[i];qHrch[i] = q1[i1]*mesh.V()[i1];}
for (i=0;i<cellsHdrn.size();i++) {i1=cellsHdrn[i];qHdrn[i]=(hSu[i1]-h()[i1]*hSp[i1])*mesh.V()[i1];}
for (i=0;i<cellsHghb.size();i++) {i1=cellsHghb[i];qHghb[i]=(hSu[i1]-h()[i1]*hSp[i1])*mesh.V()[i1];}

//phiw =  hEqn.flux() ;
//Info<<"phi "<<phiw<<endl;
//volScalarField q1 = fvc::div(phiw);
//volScalarField d1 = stor*fvc::ddt(h);
//for (ilay=0;ilay<nlay;ilay++) {
	//i1=ilay*ncell_lay+ic0;
	//Info<<"h "<<h()[i1]<<" v "<<mesh.V()[i1]<<" sw "<<sw[i1]<<" q1 "<<q1[i1]<<" d1 "<<d1[i1]<<" SuSp "<<hSu[i1]-h()[i1]*hSp[i1]<<endl;
	//}
